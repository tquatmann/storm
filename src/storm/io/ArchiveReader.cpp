#include "storm/io/ArchiveReader.h"

namespace storm::io {

#ifdef STORM_HAVE_LIBARCHIVE
/*!
 * Checks the result of an archive operation and throws an exception if the result is not ok.
 */
void checkResult(archive* arch, auto resultCode) {
    static_assert(ARCHIVE_OK == 0, "Expected that return value >= 0 means a valid result");
    STORM_LOG_THROW(arch != nullptr, storm::exceptions::FileIoException, "Unexpected result: Archive not loaded.");
    STORM_LOG_WARN_COND(std::cmp_greater_equal(resultCode, ARCHIVE_OK), "Unexpected result from archive: " << archive_error_string(arch) << ".");
    STORM_LOG_THROW(std::cmp_greater_equal(resultCode, ARCHIVE_WARN), storm::exceptions::FileIoException,
                    "Unexpected result from archive: " << archive_error_string(arch) << ".");
}

void ArchiveReader::ArchiveDeleter::operator()(archive* arch) const noexcept {
    if (arch) {
        // archives created for reading OR writing can be freed by archive_free()
        archive_free(arch);
    }
}

ArchiveReader::ArchiveReadEntry::ArchiveReadEntry(archive_entry* currentEntry, archive* archive) : _currentEntry(currentEntry), _archive(archive) {
    STORM_LOG_ASSERT(_currentEntry, "No valid entry loaded.");
}
#endif

std::filesystem::path ArchiveReader::ArchiveReadEntry::name() const {
#ifdef STORM_HAVE_LIBARCHIVE
    STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");
    std::filesystem::path result;
    char const* path = archive_entry_pathname(_currentEntry);
    if (path) {
        result = path;
    }
    return result;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

bool ArchiveReader::ArchiveReadEntry::isDir() const {
#ifdef STORM_HAVE_LIBARCHIVE
    STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");
    return archive_entry_filetype(_currentEntry) == AE_IFDIR;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

template<typename T>
using Vec = typename ArchiveReader::ArchiveReadEntry::VectorType<T>;

template<typename T, std::endian Endianness>
    requires(std::is_arithmetic_v<T>)
Vec<T> ArchiveReader::ArchiveReadEntry::toVector() {
#ifdef STORM_HAVE_LIBARCHIVE
    using BucketType = decltype(std::declval<storm::storage::BitVector&>().getBucket({}));
    constexpr bool IsBitVector = std::is_same_v<T, bool>;
    using DataType = std::conditional_t<IsBitVector, BucketType, T>;  // for BitVectors, we use uint64_t as the underlying type
    constexpr bool NativeEndianness = Endianness == std::endian::native;
    STORM_LOG_THROW(_currentEntry, storm::exceptions::FileIoException, "No valid entry loaded.");

    // Prepare the vector to store the data, using given size (if available)
    Vec<T> content;
    auto entrySize = archive_entry_size(_currentEntry);
    checkResult(_archive, entrySize);
    entrySize = std::max<decltype(entrySize)>(entrySize, 0);
    if constexpr (IsBitVector) {
        content.resize(entrySize * 8);  // 8 bits in a byte
    } else {
        // For other types, we reserve the number of elements
        content.reserve(entrySize / sizeof(DataType));
    }

    [[maybe_unused]] uint64_t bucketCount = 0;  // only used for BitVector content
    // Helper function to add data to the content
    auto append = [&content, &bucketCount](std::ranges::input_range auto&& data) {
        if constexpr (IsBitVector) {
            content.grow(bucketCount + data.size() * sizeof(BucketType) * 8);  // 8 bits in a byte
            for (auto bits : data) {
                content.setBucket(
                    bucketCount,
                    storm::utility::reverseBits(
                        bits));  // Our bit vectors store the items in reverse order, i.e., the first item is indicated by the most significant bit
                ++bucketCount;
            }
        } else {
            (void)bucketCount;  // silences unused lambda capture warning
            content.insert(content.end(), data.begin(), data.end());
        }
    };

    // todo: try out std::array as buffer instead
    static_assert(BufferSize % sizeof(DataType) == 0, "Buffer size should be a multiple of sizeof(DataType).");
    std::vector<char> buffer(BufferSize);
    auto bytesRead = archive_read_data(_archive, buffer.data(), BufferSize);
    checkResult(_archive, bytesRead);
    while (true) {
        // process the current buffer contents
        auto const numValues = bytesRead / sizeof(DataType);  // number of values that we can now append
        if constexpr (NativeEndianness || sizeof(DataType) == 1) {
            append(std::span<const DataType>(reinterpret_cast<const DataType*>(buffer.data()), numValues));
        } else {
            append(std::span<const DataType>(reinterpret_cast<const DataType*>(buffer.data()), numValues) |
                   std::ranges::views::transform(storm::utility::byteSwap<DataType>));
        }

        // put the next chunk into the buffer
        auto offsetBytes = bytesRead % sizeof(DataType);  // number of bytes that could not be processed in this round
        if (offsetBytes > 0 && numValues > 0) {
            // if some of the bytes could not be processed, we copy them to the beginning of the buffer for the next read
            std::copy(buffer.data() + bytesRead - offsetBytes, buffer.data() + bytesRead, buffer.data());
        }
        bytesRead = archive_read_data(_archive, buffer.data() + offsetBytes, BufferSize - offsetBytes);
        checkResult(_archive, bytesRead);
        if (bytesRead == 0) {
            STORM_LOG_THROW(
                offsetBytes == 0, storm::exceptions::FileIoException,
                "Archive entry could not be extracted as vector of a " << sizeof(DataType) << "-bytes type: " << offsetBytes << " bytes left in the buffer.");
            break;  // no more data to read
        }
        bytesRead += offsetBytes;  // actual number of bytes to process in the buffer
    }

    // Resize the content to the actual size
    if constexpr (IsBitVector) {
        content.resize(bucketCount * sizeof(BucketType) * 8);  // 8 bits in a byte
    } else {
        content.shrink_to_fit();
    }

    // We have read the data, i.e. it is no longer readable
    _archive = nullptr;
    return content;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

std::string ArchiveReader::ArchiveReadEntry::toString() {
#ifdef STORM_HAVE_LIBARCHIVE
    // Prepare the vector to store the data, using given size (if available)
    std::string content;
    auto const entrySize = archive_entry_size(_currentEntry);
    checkResult(_archive, entrySize);
    content.reserve(std::max<decltype(entrySize)>(entrySize, 0));

    // todo: try out std::array as buffer instead
    std::vector<char> buffer(BufferSize);
    la_ssize_t bytesRead = 0;
    while ((bytesRead = archive_read_data(_archive, buffer.data(), BufferSize)) > 0) {
        content.append(buffer.data(), bytesRead);
    }
    STORM_LOG_THROW(bytesRead >= 0, storm::exceptions::FileIoException, "Failed to read data from archive. " << archive_error_string(_archive) << ".");
    content.shrink_to_fit();

    // We have read the data, i.e. it is no longer readable
    _archive = nullptr;
    return content;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

#ifdef STORM_HAVE_LIBARCHIVE
ArchiveReader::Iterator::Iterator(std::filesystem::path const& filename) : _archive(archive_read_new(), ArchiveDeleter{}), _currentEntry(nullptr) {
    STORM_LOG_THROW(_archive, storm::exceptions::FileIoException, "Failed to create archive reader.");
    // Enable all filters (e.g., gzip, bzip2, xz) and all formats (tar, zip, etc.)
    checkResult(_archive.get(), archive_read_support_filter_all(_archive.get()));
    checkResult(_archive.get(), archive_read_support_format_all(_archive.get()));
    // A typical block size of 10240 is recommended by libarchive documentation
    checkResult(_archive.get(), archive_read_open_filename(_archive.get(), filename.c_str(), 10240));
    ++*this;  // Move to the first entry
}
#endif

bool ArchiveReader::Iterator::operator==(Iterator const& other) const {
#ifdef STORM_HAVE_LIBARCHIVE
    return _currentEntry == other._currentEntry;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

bool ArchiveReader::Iterator::operator!=(Iterator const& other) const {
#ifdef STORM_HAVE_LIBARCHIVE
    return _currentEntry != other._currentEntry;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

/*!
 * Move to the next entry in the archive.
 */
typename ArchiveReader::Iterator& ArchiveReader::Iterator::operator++() {
#ifdef STORM_HAVE_LIBARCHIVE
    int r = archive_read_next_header(_archive.get(), &_currentEntry);
    if (r == ARCHIVE_EOF) {
        // End of archive
        _currentEntry = nullptr;
        _archive.reset();
    } else {
        checkResult(_archive.get(), r);
    }
    return *this;
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

typename ArchiveReader::ArchiveReadEntry ArchiveReader::Iterator::operator*() const {
#ifdef STORM_HAVE_LIBARCHIVE
    return ArchiveReadEntry(_currentEntry, _archive.get());
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

ArchiveReader::ArchiveReader(std::filesystem::path const& file) : file(file) {};

typename ArchiveReader::Iterator ArchiveReader::begin() const {
#ifdef STORM_HAVE_LIBARCHIVE
    return Iterator(file);
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

typename ArchiveReader::Iterator ArchiveReader::end() const {
#ifdef STORM_HAVE_LIBARCHIVE
    return Iterator();
#else
    throw storm::exceptions::NotSupportedException() << "Reading archives is not supported. Storm is compiled without LibArchive.";
#endif
}

ArchiveReader openArchive(std::filesystem::path const& file) {
    return ArchiveReader(file);
}

template Vec<char> ArchiveReader::ArchiveReadEntry::toVector<char, std::endian::little>();
template Vec<bool> ArchiveReader::ArchiveReadEntry::toVector<bool, std::endian::little>();
template Vec<uint32_t> ArchiveReader::ArchiveReadEntry::toVector<uint32_t, std::endian::little>();
template Vec<uint64_t> ArchiveReader::ArchiveReadEntry::toVector<uint64_t, std::endian::little>();
template Vec<int64_t> ArchiveReader::ArchiveReadEntry::toVector<int64_t, std::endian::little>();
template Vec<double> ArchiveReader::ArchiveReadEntry::toVector<double, std::endian::little>();

}  // namespace storm::io