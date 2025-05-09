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
    void addFile(std::string const& archivePath, char const* data, std::size_t size) {
        archive_entry* entry = archive_entry_new();
        STORM_LOG_THROW(entry, storm::exceptions::FileIoException, "Failed to create archive entry.");
        STORM_LOG_THROW(data == nullptr || size > 0, storm::exceptions::FileIoException, "Data buffer is empty.");

        // Fill in metadata: path, file size, file type, permissions, etc.
        archive_entry_set_pathname(entry, archivePath.c_str());
        archive_entry_set_size(entry, size);
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0777);

        // Write the header (metadata) to the archive
        checkResult(archive_write_header(_archive.get(), entry), entry);

        // Write the file contents
        if (data != nullptr) {
            checkResult(archive_write_data(_archive.get(), data, size));
        }

        // Free the entry metadata after we finish writing
        archive_entry_free(entry);
    }

    template<typename T, std::endian Endianness = std::endian::little>
        requires(std::is_arithmetic_v<T> && !std::is_same_v<T, bool>)
    void addBinaryFile(std::string const& archivePath, std::vector<T> const& data) {
        auto const numBytes = data.size() * sizeof(T);
        if constexpr (Endianness == std::endian::native || sizeof(T) == 1) {
            // can write directly without a buffer
            addFile(archivePath, reinterpret_cast<char const*>(data.data()), numBytes);
        } else {
            // need to swap bytes. This uses the buffer
            _buffer.resize(std::max(_buffer.size(), numBytes));
            std::span<T> buffer_T(reinterpret_cast<T*>(_buffer.data()), data.size());
            for (std::size_t i = 0; i < data.size(); ++i) {
                buffer_T[i] = storm::utility::byteSwap(data[i]);
            }
            addFile(archivePath, _buffer.data(), numBytes);
        }
    }

    void addBinaryFile(std::string const& archivePath, storm::storage::BitVector const& data) {
        using BucketType = decltype(std::declval<storm::storage::BitVector&>().getBucket({}));
        auto const numBytes = data.bucketCount() * sizeof(BucketType);
        // need to reverse bits and potentially swap bytes
        _buffer.resize(std::max(_buffer.size(), numBytes));
        std::span<BucketType> buffer_BT(reinterpret_cast<BucketType*>(_buffer.data()), data.bucketCount());
        for (uint64_t i = 0; i < data.bucketCount(); ++i) {
            // reverse bits so that the first bit of the first byte is data.get(0).
            buffer_BT[i] = storm::utility::reverseBits<BucketType, std::endian::native == std::endian::little>(data.getBucket(i));
        }
        addFile(archivePath, _buffer.data(), numBytes);
    }

    void addTextFile(std::string const& archivePath, std::string const& data) {
        addFile(archivePath, data.c_str(), data.size());
    }

   private:
    /*!
     * Checks the result of an archive operation and throws an exception if the result is not ok.
     */
    void checkResult(int resultCode, archive_entry* entry = nullptr) {
        STORM_LOG_WARN_COND(resultCode >= ARCHIVE_OK, "Unexpected result from archive: " << archive_error_string(_archive.get()) << ".");
        if (resultCode < ARCHIVE_WARN) {
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

    std::unique_ptr<archive, ArchiveDeleter> _archive;
    std::vector<char> _buffer;
};
}  // namespace storm::io
