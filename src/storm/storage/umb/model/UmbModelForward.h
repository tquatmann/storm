#pragma once

namespace storm::umb {

class UmbModelBase;

enum class StorageType;

template<StorageType Storage>
class UmbModel;

struct ModelIndex;

}  // namespace storm::umb