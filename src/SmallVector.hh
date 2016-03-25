#pragma once
#include "SmallMatrix.hh"
#include "VectorOf.hh"

namespace linalgwrap {
/** Using statement to define a SmallVector */
template <typename Scalar>
using SmallVector = VectorOf<SmallMatrix<Scalar>>;
}  // end namespace linalgwrap
