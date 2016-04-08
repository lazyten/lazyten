#pragma once
#include "ArmadilloMatrix.hh"

namespace linalgwrap {
#if defined LINALGWRAP_HAVE_ARMADILLO
// Forward declare:
template <typename Scalar>
class ArmadilloMatrix;

// Use armadillo if available
template <typename Scalar>
using SmallMatrix = ArmadilloMatrix<Scalar>;
#else
template <typename Scalar>
class SmallMatrix {
    static_assert(false, "No default implementation for SmallMatrix.");
};
#endif

}  // end namespace linalgwrap
