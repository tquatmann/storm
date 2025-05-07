#pragma once

#include <archive.h>
#include <archive_entry.h>

#include <bit>
#include <concepts>
#include <memory>
#include <string>
#include <vector>

#include "storm/exceptions/FileIoException.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace detail {

/*!
 * Auxiliary struct to delete archive objects (not the archive from disk!)
 */
struct ArchiveDeleter {
    void operator()(archive* arch) const noexcept {
        if (arch) {
            // archives created for reading OR writing can be freed by archive_free()
            archive_free(arch);
        }
    }
};

class ArchiveEntryReadIterator {
   public:
    ArchiveEntryReadIterator() = default;

    ArchiveEntryReadIterator(std::string const& filename) : _archive(archive_read_new(), detail::ArchiveDeleter{}), _currentEntry(nullptr) {
        STORM_LOG_THROW(_archive, storm::exceptions::FileIoException, "Failed to create archive reader.");
        // Enable all filters (e.g., gzip, bzip2, xz) and all formats (tar, zip, etc.)
        checkResult(archive_read_support_filter_all(_archive.get()));
        checkResult(archive_read_support_format_all(_archive.get()));
        // A typical block size of 10240 is recommended by libarchive documentation
        checkResult(archive_read_open_filename(_archive.get(), filename.c_str(), 10240));
    }

    bool operator==(ArchiveEntryReadIterator const& other) const = default;

    /*!
     * Move to the next entry in the archive. Frees the archive
     */
    ArchiveEntryReadIterator& operator++() {
        int r = archive_read_next_header(_archive.get(), &_currentEntry);
        if (r == ARCHIVE_EOF) {
            // End of archive
            _currentEntry = nullptr;
            _archive.reset();
        }
        checkResult(r);
        return *this;
    }

    /*!
     * Get the current entry’s path (filename) inside the archive.
     * Returns an empty string if the entry does not exist.
     */
    std::string name() const {
        STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");
        const char* path = archive_entry_pathname(_currentEntry);
        return path ? path : std::string{};
    }

    template<typename T>
    using VectorType = std::conditional_t<std::is_same_v<T, bool>, storm::storage::BitVector, std::vector<T>>;

    /*!
     * extracts the current entry’s data as a vector of the given type.
     */
    template<typename T, std::endian Endianness = std::endian::little>
        requires(std::is_arithmetic_v<T>)
    VectorType<T> dataAsVector() {
        constexpr bool NativeEndianness = Endianness == std::endian::native;
        constexpr auto BytesPerEntry = std::is_same_v<T, bool> ? sizeof(uint64_t) : sizeof(T);
        STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");

        // Prepare the vector to store the data, using given size (if available)
        VectorType<T> content;
        [[maybe_unused]] uint64_t bucketCount = 0;  // only used for BitVector content
        std::cout << "Reading entry: " << name() << " (size=" << archive_entry_size(_currentEntry) << " bytes)\n";
        if constexpr (std::is_same_v<T, bool>) {
            content.resize(archive_entry_size(_currentEntry) * 8);  // 8 bits in a byte
        } else {
            // For other types, we reserve the number of elements
            content.reserve(archive_entry_size(_currentEntry) / BytesPerEntry);
        }

        // Helper function to add data to the content
        auto append = [&content, &bucketCount](std::ranges::input_range auto&& data) {
            if constexpr (std::is_same_v<T, bool>) {
                content.grow(bucketCount + data.size() * 64);  // 64 bits in a bucket
                for (auto bits : data) {
                    content.setBucket(
                        bucketCount,
                        storm::utility::reverseBits(
                            bits));  // Our bit vectors store the items in reverse order, i.e., the first item is indicated by the most significant bit
                    ++bucketCount;
                }
            } else {
                content.insert(content.end(), data.begin(), data.end());
            }
        };

        // todo: try out std::array as buffer instead
        constexpr size_t bufferSize = 8192;
        std::vector<char> buffer(bufferSize);
        la_ssize_t bytesRead = 0;
        la_ssize_t offsetBytes = 0;  // used in case the number of bytes read is not a multiple of BytesPerEntry
        while ((bytesRead = archive_read_data(_archive.get(), buffer.data() + offsetBytes, bufferSize - offsetBytes)) > 0) {
            bytesRead += offsetBytes;                          // actual number of bytes to process in the buffer
            offsetBytes = bytesRead % BytesPerEntry;           // number of bytes that can not be processed in this round
            auto const numValues = bytesRead / BytesPerEntry;  // values that we can now append
            if constexpr (NativeEndianness || BytesPerEntry == 1) {
                append(std::span<const T>(reinterpret_cast<const T*>(buffer.data()), numValues));
            } else {
                append(std::span<const T>(reinterpret_cast<const T*>(buffer.data()), numValues) | std::ranges::views::transform(storm::utility::byteSwap<T>));
            }
            if (offsetBytes > 0 && numValues > 0) {
                std::cout << " have offset bytes: " << offsetBytes << " bytes\n";
                // if some of the bytes could not be processed, we copy them to the beginning of the buffer for the next read
                std::copy(buffer.data() + bytesRead - offsetBytes, buffer.data() + bytesRead, buffer.data());
            }
        }
        STORM_LOG_THROW(bytesRead >= 0, storm::exceptions::FileIoException,
                        "Failed to read data from archive. " << archive_error_string(_archive.get()) << ".");
        STORM_LOG_THROW(
            offsetBytes == 0, storm::exceptions::FileIoException,
            "Archive entry could not be extracted as vector of a " << BytesPerEntry << "-bytes type: " << offsetBytes << " bytes left in the buffer.");
        std::cout << "\t actual size is " << (content.size() * BytesPerEntry) << " bytes.\n";

        // Resize the content to the actual size
        if constexpr (std::is_same_v<T, bool>) {
            content.resize(bucketCount * 64);  // 64 bits in a bucket
        } else {
            content.shrink_to_fit();
        }
        return content;
    }

   private:
    void checkResult(int resultCode) {
        STORM_LOG_WARN_COND(resultCode >= ARCHIVE_OK, "Unexpected result from archive: " << archive_error_string(_archive.get()) << ".");
        STORM_LOG_THROW(resultCode >= ARCHIVE_WARN, storm::exceptions::FileIoException,
                        "Unexpected result from  archive: " << archive_error_string(_archive.get()) << ".");
    }

   private:
    archive_entry* _currentEntry;
    std::unique_ptr<archive, detail::ArchiveDeleter> _archive;
};

class ArchiveReadHandle {
   public:
    using Iterator = ArchiveEntryReadIterator;

    ArchiveReadHandle(std::string const& filename) : filename(filename){};

    Iterator begin() const {
        return Iterator(filename);
    }

    Iterator end() const {
        return Iterator();
    }

   private:
    std::string const filename;
};

}  // namespace detail

/*!
 * Reads an archive file
 * @param filename
 * @return A range-like object to iterate over the entries in the archive.
 */
detail::ArchiveReadHandle openArchive(std::string const& filename) {
    return detail::ArchiveReadHandle(filename);
}

// ---------------------------------------------------------------------------------
// Writer wrapper class
// ---------------------------------------------------------------------------------
class ArchiveWriter {
   public:
    /// Create a new archive and open it as a file on disk. By default, it uses
    /// pax-restricted (POSIX) format. You could call other setter methods for
    /// format if needed.
    explicit ArchiveWriter(const std::string& filename) : m_archive(archive_write_new(), ArchiveDeleter{}) {
        if (!m_archive) {
            throw std::bad_alloc();
        }
        // For example, set format to TAR with restricted pax extensions
        checkResult(archive_write_set_format_pax_restricted(m_archive.get()));

        // Open file for writing
        checkResult(archive_write_open_filename(m_archive.get(), filename.c_str()));
    }

    /// Add a file to the archive. The file contents come from the data buffer.
    /// The file’s path inside the archive is given by "archivePath".
    void addFile(const std::string& archivePath, const std::vector<char>& data) {
        // Create a new entry metadata object
        archive_entry* entry = archive_entry_new();
        if (!entry) {
            throw std::bad_alloc();
        }

        // Fill in metadata: path, file size, file type, permissions, etc.
        archive_entry_set_pathname(entry, archivePath.c_str());
        archive_entry_set_size(entry, data.size());
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0644);

        // Write the header (metadata) to the archive
        int r = archive_write_header(m_archive.get(), entry);
        if (r != ARCHIVE_OK) {
            archive_entry_free(entry);
            throw ArchiveException(archive_error_string(m_archive.get()));
        }

        // Write the file contents
        if (!data.empty()) {
            la_ssize_t writtenSize = archive_write_data(m_archive.get(), data.data(), data.size());
            if (writtenSize < 0) {
                archive_entry_free(entry);
                throw ArchiveException(archive_error_string(m_archive.get()));
            }
        }

        // Free the entry metadata after we finish writing
        archive_entry_free(entry);
    }

   private:
    void checkResult(int resultCode) {
        if (resultCode != ARCHIVE_OK && resultCode != ARCHIVE_WARN) {
            throw ArchiveException(archive_error_string(m_archive.get()));
        }
    }

   private:
    std::unique_ptr<archive, ArchiveDeleter> m_archive;
};

// ---------------------------------------------------------------------------------
// Example usage in a main function
// ---------------------------------------------------------------------------------
int main() {
    try {
        // 1. Reading an archive
        ArchiveReader reader;
        reader.openFile("test.tar.gz");
        while (reader.nextEntry()) {
            std::cout << "Found entry: " << reader.currentEntryName() << "\n";
            auto data = reader.readData();
            std::cout << "  -> Entry size: " << data.size() << " bytes.\n";
        }

        // 2. Creating a new archive and adding a file
        ArchiveWriter writer("example_out.tar");
        const std::vector<char> fileContent = {'H', 'e', 'l', 'l', 'o', ' ', 'f', 'r', 'o', 'm', ' ', 'C', '+', '+', '\n'};
        writer.addFile("hello.txt", fileContent);

        std::cout << "Archive created successfully!\n";
    } catch (const ArchiveException& e) {
        std::cerr << "Archive error: " << e.what() << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Standard exception: " << e.what() << '\n';
    }
    return 0;
}