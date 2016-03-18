#pragma once
#include <type_traits>
#include "Matrix_i.hh"
#include "StoredMatrix_i.hh"

namespace linalgwrap {

template <typename T>
class IsMatrix : public std::is_base_of<Matrix_i<typename T::scalar_type>, T> {
};

template <typename Matrix>
class IsStoredMatrix
      : public std::is_base_of<StoredMatrix_i<typename Matrix::scalar_type>,
                               Matrix> {};
}  // namespace linalgwrap
