#pragma once
#include <filesystem>
#include <ranges>
#include <span>

#include <boost/core/typeinfo.hpp>
#include <boost/iostreams/device/mapped_file.hpp>

#include "storm/exceptions/FileIoException.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"

namespace storm::io {
static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big, "This code is not supported for mixed endian systems.");

template<typename T>
concept IsBinaryFileViewable = std::is_arithmetic_v<T> && !std::is_same_v<T, bool>;

template<typename T, std::endian Endianness = std::endian::little>
    requires IsBinaryFileViewable<T>
class BinaryFileViewer {
   public:
    static const bool NativeEndianness = Endianness == std::endian::native;
    using value_type = T;

    BinaryFileViewer(std::filesystem::path const& filename) : file(createMappedFile(filename)), rangeView(createRangeView(file)) {}

    const T& operator[](std::size_t index) const {
        return rangeView[index];
    }

    std::size_t size() const {
        return rangeView.size();
    }

    auto begin() const {
        return rangeView.begin();
    }

    auto end() const {
        return rangeView.end();
    }

    auto data() const
        requires(NativeEndianness || sizeof(T) == 1)
    {
        return rangeView.data();
    }

   private:
    using MappedFileType = boost::iostreams::mapped_file_source;
    static MappedFileType createMappedFile(std::filesystem::path const& filename) {
        if (std::filesystem::is_empty(filename)) {
            return MappedFileType();
        }
        MappedFileType file(filename);
        STORM_LOG_THROW(file.is_open(), storm::exceptions::FileIoException, "Could not open file '" << filename << "'.");
        STORM_LOG_THROW(file.size() % sizeof(T) == 0, storm::exceptions::FileIoException,
                        "Error when opening file '" << filename << "'. File size " << file.size() << " is not a multiple of " << sizeof(T) << ".");
        STORM_LOG_DEBUG("Opened file '" << filename << "' (size=" << file.size() << " bytes, type='" << boost::core::demangled_name(BOOST_CORE_TYPEID(T))
                                        << ", length=" << (file.size() / sizeof(T)) << ").");
        return file;
    }

    static auto createRangeView(MappedFileType const& file) {
        const T* data = file.is_open() ? reinterpret_cast<const T*>(file.data()) : nullptr;
        std::size_t size = file.is_open() ? file.size() / sizeof(T) : 0;
        if constexpr (NativeEndianness || sizeof(T) == 1) {
            return std::span<const T>(data, size);
        } else {
            return std::span<const T>(data, size) | std::ranges::views::transform(storm::utility::byteSwap<T>);
        }
    }

    using RangeViewType = std::invoke_result_t<decltype(createRangeView), MappedFileType const&>;
    boost::iostreams::mapped_file_source file;
    RangeViewType rangeView;
};
static_assert(std::ranges::random_access_range<BinaryFileViewer<double>>);
static_assert(std::ranges::random_access_range<BinaryFileViewer<uint64_t>>);
static_assert(std::ranges::random_access_range<BinaryFileViewer<int32_t>>);

}  // namespace storm::io