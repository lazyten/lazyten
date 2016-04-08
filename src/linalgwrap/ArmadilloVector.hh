#pragma once
#include "linalgwrap/ArmadilloMatrix.hh"
#include "linalgwrap/VectorOf.hh"

namespace linalgwrap {
#ifdef LINALGWRAP_HAVE_ARMADILLO
template <typename Scalar>
using ArmadilloVector = VectorOf<ArmadilloMatrix<Scalar>>;
#endif
}  // end namespace linalgwrap
