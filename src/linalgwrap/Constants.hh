#ifndef LIBLINALG_CONSTANTS_HH_
#define LIBLINALG_CONSTANTS_HH_
#include <limits>

namespace linalgwrap {

template <typename Scalar>
struct Constants {
    typedef Scalar scalar_type;

    //! Constant representing 0
    static constexpr scalar_type zero = scalar_type(0.);

    //! Constant representing 1
    static constexpr scalar_type one = scalar_type(1.);

    /** \brief Constant which should be used for invalid values
     *
     * \note This value may well be a valid ``scalar_type`` value, so
     * do not check against it in order to see if the ``scalar_type`` value
     * is invalid.
     */
    static constexpr scalar_type invalid =
          std::numeric_limits<scalar_type>::has_signaling_NaN
                ? std::numeric_limits<scalar_type>::signaling_NaN()
                : std::numeric_limits<scalar_type>::max();
};
}

#endif
