#pragma once
#include <ArmadilloMatrix.hh>

namespace linalgwrap {
#if defined LINALGWRAP_HAVE_ARMADILLO
// Use armadillo if available
template <typename Scalar>
using SmallMatrix = ArmadilloMatrix<Scalar>;
#else
static_assert(false, "No implementation for SmallMatrix.");
#endif
}  // end namespace linalgwrap
