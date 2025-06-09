#pragma once

#include <archive.h>
#include <archive_entry.h>

#include <bit>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "storm/exceptions/FileIoException.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::io {
static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big, "This code is not supported for mixed endian systems.");

template<typename R>
concept WritableRange = std::ranges::input_range<R> && std::is_arithmetic_v<std::ranges::range_value_t<R>>;

template<typename R, std::endian Endianness>
concept WritableWithoutBuffer =
    WritableRange<R> && std::same_as<std::remove_cvref_t<R>, std::vector<std::ranges::range_value_t<R>>> &&
    !std::is_same_v<std::ranges::range_value_t<R>, bool> && (Endianness == std::endian::native || sizeof(std::ranges::range_value_t<R>) == 1);

/*!
  Create a new archive and open it as a file on disk.
 */
class ArchiveWriter {
   public:
    explicit ArchiveWriter(std::string const& filename) : _archive(archive_write_new(), ArchiveDeleter{}) {
        STORM_LOG_THROW(_archive, storm::exceptions::FileIoException, "Failed to create archive reader.");

        // Set format to gzipped TAR with restricted pax extensions
        archive_write_add_filter_gzip(_archive.get());
        checkResult(archive_write_set_format_pax_restricted(_archive.get()));

        // Open file for writing
        checkResult(archive_write_open_filename(_archive.get(), filename.c_str()));
    }

    /*!
      Adds a (sub-) directory to the archive
     */
    void addDirectory(std::string const& archivePath) {
        archive_entry* entry = archive_entry_new();
        STORM_LOG_THROW(entry, storm::exceptions::FileIoException, "Failed to create archive entry.");

        // Fill in metadata: path, file type, permissions, etc.
        archive_entry_set_pathname(entry, archivePath.c_str());
        archive_entry_set_filetype(entry, AE_IFDIR);
        archive_entry_set_perm(entry, 0777);

        // Write the header (metadata) to the archive
        checkResult(archive_write_header(_archive.get(), entry), entry);

        // Free the entry metadata after we finish writing
        archive_entry_free(entry);
    }

    /*!
     * Add a file to the archive. The file contents come from the data buffer.
     * The file’s path inside the archive is given by "archivePath".
     */
    void addFile(std::string const& archivePath, char const* data, std::size_t const size) {
        std::size_t startOfChunk = 0;
        auto getNextChunk = [&]() {
            auto const chunkSize = std::min<std::size_t>(size - startOfChunk, BufferSize);
            auto const chunk = std::span<char const>(data + startOfChunk, chunkSize);
            startOfChunk += chunkSize;
            return chunk;
        };
        addFileFromChunks(archivePath, getNextChunk, size);
    }

    template<WritableRange Range, std::endian Endianness = std::endian::little>
        requires(WritableWithoutBuffer<Range, Endianness>)
    void addBinaryFile(std::string const& archivePath, Range&& data) {
        auto const numBytes = data.size() * sizeof(std::ranges::range_value_t<Range>);
        addFile(archivePath, reinterpret_cast<char const*>(data.data()), numBytes);
    }

    template<WritableRange Range, std::endian Endianness = std::endian::little>
        requires(!WritableWithoutBuffer<Range, Endianness> && !std::is_same_v<std::ranges::range_value_t<Range>, bool>)
    void addBinaryFile(std::string const& archivePath, Range const& data) {
        using T = std::ranges::range_value_t<Range>;
        auto const numBytes = data.size() * sizeof(T);
        static_assert(BufferSize % sizeof(T) == 0, "Buffer size must be a multiple of sizeof(T).");
        std::vector<T> buffer(BufferSize / sizeof(T));  // todo check array

        uint64_t startOfChunk = 0;
        auto getNextChunk = [&data, &buffer, &startOfChunk]() {
            auto const endOfChunk = std::min<uint64_t>(startOfChunk + buffer.size(), data.size());
            auto bufferIt = buffer.begin();
            for (uint64_t i = startOfChunk; i < endOfChunk; ++i) {
                if constexpr (std::endian::native != Endianness && sizeof(T) > 1) {
                    // copy into buffer with byteswap
                    *bufferIt = storm::utility::byteSwap(data[i]);
                } else {
                    // copy into buffer
                    *bufferIt = data[i];
                }
                ++bufferIt;
            }
            std::span<const char> chunk(reinterpret_cast<const char*>(buffer.data()), (endOfChunk - startOfChunk) * sizeof(T));
            startOfChunk = endOfChunk;
            return chunk;
        };
        addFileFromChunks(archivePath, getNextChunk, numBytes);
    }

    void addBinaryFile(std::string const& archivePath, storm::storage::BitVector const& data) {
        using BucketType = decltype(std::declval<storm::storage::BitVector&>().getBucket({}));
        static_assert(BufferSize % sizeof(BucketType) == 0, "Buffer size must be a multiple of sizeof(BucketType).");
        std::vector<BucketType> buffer(BufferSize / sizeof(BucketType));  // todo check array

        // need to reverse bits and potentially swap bytes
        uint64_t startOfChunk = 0;
        auto getNextChunk = [&data, &buffer, &startOfChunk]() {
            auto const endOfChunk = std::min<uint64_t>(startOfChunk + buffer.size(), data.bucketCount());
            auto bufferIt = buffer.begin();
            for (uint64_t i = startOfChunk; i < endOfChunk; ++i) {
                // reverse bits so that the first bit of the first byte is data.get(0).
                bool constexpr NativeLittleEndian = std::endian::native == std::endian::little;
                *bufferIt = storm::utility::reverseBits<BucketType, NativeLittleEndian>(data.getBucket(i));
                ++bufferIt;
            }
            std::span<const char> chunk(reinterpret_cast<const char*>(buffer.data()), (endOfChunk - startOfChunk) * sizeof(BucketType));
            startOfChunk = endOfChunk;
            return chunk;
        };
        uint64_t const numBytes = data.bucketCount() * sizeof(BucketType);
        addFileFromChunks(archivePath, getNextChunk, numBytes);
    }

    void addTextFile(std::string const& archivePath, std::string const& data) {
        addFile(archivePath, data.c_str(), data.size());
    }

   private:
    /*!
     * Checks the result of an archive operation and throws an exception if the result is not ok.
     */
    void checkResult(auto resultCode, archive_entry* entry = nullptr) {
        STORM_LOG_THROW(_archive, storm::exceptions::FileIoException, "Unexpected result: Archive not loaded.");

        STORM_LOG_WARN_COND(std::cmp_greater_equal(resultCode, ARCHIVE_OK), "Unexpected result from archive: " << archive_error_string(_archive.get()) << ".");
        if (std::cmp_less(resultCode, ARCHIVE_WARN)) {
            if (entry) {
                archive_entry_free(entry);
            }
            STORM_LOG_THROW(false, storm::exceptions::FileIoException, "Unexpected result from archive: " << archive_error_string(_archive.get()) << ".");
        }
    }

   private:
    struct ArchiveDeleter {
        void operator()(archive* arch) const noexcept {
            if (arch) {
                // archives created for reading OR writing can be freed by archive_free()
                archive_free(arch);
            }
        }
    };

    /*!
     * Add a file to the archive.
     * @param archivePath the path inside the archive
     * @param getNextChunk returns the file contents in bytes, divided into multiple chunks
     * @param size the size of the file in bytes
     */
    template<typename F>
    void addFileFromChunks(std::string const& archivePath, F getNextChunk, size_t const size) {
        archive_entry* entry = archive_entry_new();
        STORM_LOG_THROW(entry, storm::exceptions::FileIoException, "Failed to create archive entry.");
        STORM_LOG_THROW(size > 0, storm::exceptions::FileIoException, "Data buffer is empty.");

        // Fill in metadata: path, file size, file type, permissions, etc.
        archive_entry_set_pathname(entry, archivePath.c_str());
        archive_entry_set_size(entry, size);
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0777);

        // Write the header (metadata) to the archive
        checkResult(archive_write_header(_archive.get(), entry), entry);

        // Write the file contents
        uint64_t bytesWritten = 0;
        for (auto chunk = getNextChunk(); chunk.size() > 0 && bytesWritten < size; chunk = getNextChunk()) {
            auto res = archive_write_data(_archive.get(), chunk.data(), chunk.size());
            checkResult(res, entry);
            bytesWritten += res;
        }
        STORM_LOG_WARN_COND(bytesWritten == size, "When writing file '" << archivePath << "' to archive, " << bytesWritten << " bytes were written but " << size
                                                                        << " bytes were expected.");

        // Free the entry metadata after we finish writing
        archive_entry_free(entry);
    }

    std::unique_ptr<archive, ArchiveDeleter> _archive;
    static constexpr size_t BufferSize = 8192;
};
}  // namespace storm::io
