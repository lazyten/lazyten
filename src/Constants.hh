#ifndef LIBLINALG_CONSTANTS_HH_
#define LIBLINALG_CONSTANTS_HH_

namespace linalgwrap {

template <typename Scalar>
struct Constants {
    typedef Scalar scalar_type;

    static constexpr scalar_type zero = scalar_type(0.);
    static constexpr scalar_type one = scalar_type(1.);
};
}

#endif
