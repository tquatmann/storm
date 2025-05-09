#pragma once

#include <archive.h>
#include <archive_entry.h>

#include <bit>
#include <memory>
#include <string>
#include <vector>

#include "storm/exceptions/FileIoException.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::io {
namespace internal {

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

/*!
 * Checks the result of an archive operation and throws an exception if the result is not ok.
 */
void checkResult(archive* arch, int resultCode) {
    STORM_LOG_WARN_COND(resultCode >= ARCHIVE_OK, "Unexpected result from archive: " << archive_error_string(arch) << ".");
    STORM_LOG_THROW(resultCode >= ARCHIVE_WARN, storm::exceptions::FileIoException, "Unexpected result from archive: " << archive_error_string(arch) << ".");
}

/*!
 * Object that reads the archive entry.
 */
class ArchiveReadEntry {
   public:
    ArchiveReadEntry(archive_entry* currentEntry, archive* archive) : _currentEntry(currentEntry), _archive(archive) {
        STORM_LOG_ASSERT(_currentEntry, "No valid entry loaded.");
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

    bool isDir() const {
        STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");
        return archive_entry_filetype(_currentEntry) == AE_IFDIR;
    }

    bool isReadable() const {
        return _archive == nullptr;
    }

    template<typename T>
    using VectorType = std::conditional_t<std::is_same_v<T, bool>, storm::storage::BitVector, std::vector<T>>;

    /*!
     * Extracts the current entry’s data as a vector of the given type.
     * @note this consumes the archive entry contents, i.e., calling this can only be called once
     */
    template<typename T, std::endian Endianness = std::endian::little>
        requires(std::is_arithmetic_v<T>)
    VectorType<T> toVector() {
        using BucketType = decltype(std::declval<storm::storage::BitVector&>().getBucket({}));
        using DataType = std::conditional_t<std::is_same_v<T, bool>, BucketType, T>;  // for BitVectors, we use uint64_t as the underlying type
        constexpr bool NativeEndianness = Endianness == std::endian::native;
        STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");

        // Prepare the vector to store the data, using given size (if available)
        VectorType<T> content;
        [[maybe_unused]] uint64_t bucketCount = 0;  // only used for BitVector content
        std::cout << "Reading entry: " << name() << " (size=" << archive_entry_size(_currentEntry) << " bytes, " << sizeof(DataType) << " bytes per entry)\n";
        if constexpr (std::is_same_v<T, bool>) {
            content.resize(archive_entry_size(_currentEntry) * 8);  // 8 bits in a byte
        } else {
            // For other types, we reserve the number of elements
            content.reserve(archive_entry_size(_currentEntry) / sizeof(DataType));
        }

        // Helper function to add data to the content
        auto append = [&content, &bucketCount](std::ranges::input_range auto&& data) {
            if constexpr (std::is_same_v<T, bool>) {
                content.grow(bucketCount + data.size() * sizeof(BucketType) * 8);  // 8 bits in a byte
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
        while ((bytesRead = archive_read_data(_archive, buffer.data() + offsetBytes, bufferSize - offsetBytes)) > 0) {
            bytesRead += offsetBytes;                             // actual number of bytes to process in the buffer
            offsetBytes = bytesRead % sizeof(DataType);           // number of bytes that can not be processed in this round
            auto const numValues = bytesRead / sizeof(DataType);  // values that we can now append
            if constexpr (NativeEndianness || sizeof(DataType) == 1) {
                append(std::span<const DataType>(reinterpret_cast<const DataType*>(buffer.data()), numValues));
            } else {
                append(std::span<const DataType>(reinterpret_cast<const DataType*>(buffer.data()), numValues) |
                       std::ranges::views::transform(storm::utility::byteSwap<DataType>));
            }
            if (offsetBytes > 0 && numValues > 0) {
                std::cout << " have offset bytes: " << offsetBytes << " bytes\n";
                // if some of the bytes could not be processed, we copy them to the beginning of the buffer for the next read
                std::copy(buffer.data() + bytesRead - offsetBytes, buffer.data() + bytesRead, buffer.data());
            }
        }
        STORM_LOG_THROW(bytesRead >= 0, storm::exceptions::FileIoException, "Failed to read data from archive. " << archive_error_string(_archive) << ".");
        STORM_LOG_THROW(
            offsetBytes == 0, storm::exceptions::FileIoException,
            "Archive entry could not be extracted as vector of a " << sizeof(DataType) << "-bytes type: " << offsetBytes << " bytes left in the buffer.");
        std::cout << "\t actual size of vector is " << (content.size()) << " entries.\n";

        // Resize the content to the actual size
        if constexpr (std::is_same_v<T, bool>) {
            content.resize(bucketCount * sizeof(BucketType) * 8);  // 8 bits in a byte
        } else {
            content.shrink_to_fit();
        }

        // We have read the data, i.e. it is no longer readable
        _archive = nullptr;
        return content;
    }

    /*!
     * extracts the current entry’s data as a string
     * @note this consumes the archive entry contents, i.e., calling this can only be called once
     */
    std::string toString() {
        // Prepare the vector to store the data, using given size (if available)
        std::string content;
        content.reserve(archive_entry_size(_currentEntry));

        // Helper function to add data to the content
        auto append = [&content](std::ranges::input_range auto&& data) { content.insert(content.end(), data.begin(), data.end()); };

        // todo: try out std::array as buffer instead
        constexpr size_t bufferSize = 8192;
        std::vector<char> buffer(bufferSize);
        la_ssize_t bytesRead = 0;
        while ((bytesRead = archive_read_data(_archive, buffer.data(), bufferSize)) > 0) {
            content.append(buffer.data(), bytesRead);
        }
        STORM_LOG_THROW(bytesRead >= 0, storm::exceptions::FileIoException, "Failed to read data from archive. " << archive_error_string(_archive) << ".");
        content.shrink_to_fit();

        // We have read the data, i.e. it is no longer readable
        _archive = nullptr;
        return content;
    }

   private:
    archive_entry* const _currentEntry;
    archive* _archive;
};

class ArchiveReadHandle {
   public:
    class Iterator {
       public:
        Iterator() = default;

        Iterator(std::string const& filename) : _archive(archive_read_new(), internal::ArchiveDeleter{}), _currentEntry(nullptr) {
            STORM_LOG_THROW(_archive, storm::exceptions::FileIoException, "Failed to create archive reader.");
            // Enable all filters (e.g., gzip, bzip2, xz) and all formats (tar, zip, etc.)
            internal::checkResult(_archive.get(), archive_read_support_filter_all(_archive.get()));
            internal::checkResult(_archive.get(), archive_read_support_format_all(_archive.get()));
            // A typical block size of 10240 is recommended by libarchive documentation
            internal::checkResult(_archive.get(), archive_read_open_filename(_archive.get(), filename.c_str(), 10240));
            ++*this;  // Move to the first entry
        }

        bool operator==(Iterator const& other) const {
            return _currentEntry == other._currentEntry;
        }

        bool operator!=(Iterator const& other) const {
            return _currentEntry != other._currentEntry;
        }

        /*!
         * Move to the next entry in the archive. Frees the archive
         */
        Iterator& operator++() {
            int r = archive_read_next_header(_archive.get(), &_currentEntry);
            if (r == ARCHIVE_EOF) {
                // End of archive
                _currentEntry = nullptr;
                _archive.reset();
            } else {
                internal::checkResult(_archive.get(), r);
            }
            return *this;
        }

        internal::ArchiveReadEntry operator*() const {
            return internal::ArchiveReadEntry(_currentEntry, _archive.get());
        }

       private:
        archive_entry* _currentEntry{nullptr};
        std::unique_ptr<archive, ArchiveDeleter> _archive;
    };

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

}  // namespace internal

using ArchiveReadEntry = internal::ArchiveReadEntry;

/*!
 * Reads an archive file
 * @param filename
 * @return A range-like object to iterate over the entries in the archive.
 */
internal::ArchiveReadHandle openArchive(std::string const& filename) {
    return internal::ArchiveReadHandle(filename);
}
}  // namespace storm::io