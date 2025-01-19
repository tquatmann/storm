#pragma once

namespace storm::dmb {

class DmbModelBase;

enum class StorageType;

template<StorageType Storage>
class DmbModel;

struct ModelIndex;

}  // namespace storm::dmb