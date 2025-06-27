#pragma once
#include <ranges>
#include <string_view>

#include "storm/storage/umb/model/VectorType.h"
#include "storm/utility/macros.h"

namespace storm::umb {

template<std::ranges::contiguous_range InputRange>
    requires std::same_as<std::ranges::range_value_t<InputRange>, char>
auto stringVectorView(InputRange&& chars, std::ranges::contiguous_range auto&& csr) {
    STORM_LOG_ASSERT(std::ranges::size(csr) > 0, "CSR must not be empty.");
    auto const numEntries = std::ranges::size(csr) - 1;
    return std::ranges::iota_view(0ull, numEntries) |
           std::ranges::views::transform([&chars, &csr](auto i) -> std::string_view { return std::string_view(chars.data() + csr[i], csr[i + 1] - csr[i]); });
}

template<std::ranges::contiguous_range InputRange>
    requires std::same_as<std::ranges::range_value_t<InputRange>, char>
auto stringVectorViewOptionalCsr(InputRange&& chars, auto&& csr) {
    STORM_LOG_ASSERT(!csr.has_value() || std::ranges::size(csr.value()) > 0, "CSR must not be empty.");
    auto const numEntries = csr.has_value() ? std::ranges::size(csr.value()) - 1 : chars.size();
    auto const charIt = chars.begin();
    return std::ranges::iota_view(0ull, numEntries) | std::ranges::views::transform([&charIt, &csr](auto i) -> std::string_view {
               if (csr.has_value()) {
                   return std::string_view(charIt + csr.value()[i], charIt + csr.value()[i + 1]);
               } else {
                   return std::string_view(charIt + i, charIt + i + 1);  // single character
               }
           });
}

}  // namespace storm::umb