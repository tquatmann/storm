#include <bit>
#include <fstream>
#include <ranges>

#include <boost/core/typeinfo.hpp>

#include "storm/exceptions/FileIoException.h"
#include "storm/utility/bitoperations.h"
#include "storm/utility/macros.h"
namespace storm::io {

static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big, "This code is not supported for mixed endian systems.");

template<typename T>
concept IsBinaryFileWritable = std::is_arithmetic_v<T> && !std::is_same_v<T, bool>;

template<typename T, std::endian Endianness = std::endian::little>
    requires IsBinaryFileWritable<T>
class BinaryFileWriter {
   public:
    static const bool NativeEndianness = Endianness == std::endian::native;
    using value_type = T;

    BinaryFileWriter(std::string const& filename) : outputStream(filename, std::ios::out | std::ios::binary) {
        STORM_LOG_THROW(outputStream, storm::exceptions::FileIoException, "Could not open file '" << filename << "' for writing.");
    }

    void write(std::ranges::input_range auto&& r) {
        std::ranges::for_each(r, [this](T const t) { write(t); });
    }

    void write(T const t) {
        if constexpr (NativeEndianness) {
            outputStream.write(reinterpret_cast<char const*>(&t), sizeof(T));
        } else {
            T const encoded = encode(t);
            outputStream.write(reinterpret_cast<char const*>(&encoded), sizeof(T));
        }
    }
    static_assert(!std::is_same_v<T, bool>, "bool is not supported yet (unclear representation)");

   private:
    static T encode(T const t)
        requires(!NativeEndianness)
    {
        return storm::utility::byteSwap(t);
    }
    std::ofstream outputStream;
};

}  // namespace storm::io