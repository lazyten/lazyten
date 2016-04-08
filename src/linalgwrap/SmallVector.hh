#pragma once
#include "linalgwrap/SmallMatrix.hh"
#include "linalgwrap/VectorOf.hh"

namespace linalgwrap {
/** Using statement to define a SmallVector */
template <typename Scalar>
using SmallVector = VectorOf<SmallMatrix<Scalar>>;
}  // end namespace linalgwrap
