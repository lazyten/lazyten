#pragma once

namespace linalgwrap {

// Forward declarations:
template <typename IteratorCore>
class MatrixIterator;

template <typename Matrix, bool Constness>
class MatrixIteratorDefaultCore;

//! The default matrix iterator
template <typename Matrix>
using DefaultMatrixIterator =
      MatrixIterator<MatrixIteratorDefaultCore<Matrix, false>>;

//! The default matrix const iterator
template <typename Matrix>
using DefaultMatrixConstIterator =
      MatrixIterator<MatrixIteratorDefaultCore<Matrix, true>>;

}  // namespace linalgwrap

// This is intentional, in order to make sure that the cyclic dependency
// between the iterators and the Matrix_i class are correctly resolved.
#include "MatrixIterator.hh"
