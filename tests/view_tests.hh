#pragma once
#include "matrix_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

/* Namespace for default tests for views */
namespace view_tests {

/** Testing library for views
 *
 * \tparam View  The view to test
 * \tparam ViewGenArg  The argument(s) needed to generate a view or
 *                     the equivalent model. Use a std::tuple or
 *                     similar if more than one argument is required
 * */
template <typename View, typename ViewGenArg>
class TestingLibrary {
  public:
    typedef View view_type;
    // typedef typename view_type::inner_matrix_type inner_matrix_type;
    typedef typename view_type::stored_matrix_type stored_matrix_type;
    typedef typename view_type::size_type size_type;
    typedef typename view_type::scalar_type scalar_type;
    typedef ViewGenArg genarg_type;

    /** Construct a testing library instance in order to check
     *  all basic functionality of a matrix view.
     *
     *  \param view_generator   Function which takes
     *  \param model_generator
     *  \param prefix A prefix to use in the rapidcheck description string
     *  \param tolerance Numeric tolerance used for comparison in (some) checks.
     *                   It is not used for the check of basic equivalence,
     *                   where we require an agreement of 10*machine_epsilon.
     */
    TestingLibrary(
          std::function<genarg_type(void)> args_generator,
          std::function<view_type(genarg_type)> view_generator,
          std::function<stored_matrix_type(genarg_type)> model_generator,
          std::string prefix = "",
          double tolerance = TestConstants::default_num_tol);

    bool run_checks() const;

  private:
    // The testing library we use
    typedef matrix_tests::ComparativeTests<stored_matrix_type, view_type>
          comptests;

    // The caller type we use
    typedef typename comptests::template RapidcheckTestableGenerator<ViewGenArg>
          callgen_type;

    std::string m_prefix;
    callgen_type m_gen;
};

//
// ---------------------------------------
//

template <typename View, typename ViewGenArg>
TestingLibrary<View, ViewGenArg>::TestingLibrary(
      std::function<genarg_type(void)> args_generator,
      std::function<view_type(genarg_type)> view_generator,
      std::function<stored_matrix_type(genarg_type)> model_generator,
      std::string prefix, double tolerance)
      : m_prefix{prefix},
        m_gen{args_generator, view_generator, model_generator, tolerance} {}

template <typename View, typename ViewGenArg>
bool TestingLibrary<View, ViewGenArg>::run_checks() const {
    const double eps = std::numeric_limits<scalar_type>::epsilon();

    bool res = true;

    /* TODO enable if read-write views are implemented
    if (View::view_of_stored_matrix) {
        res &= rc::check(m_prefix + "Test copying stored matrix views",
                         m_gen.generate(comptests::test_copy));
    }
    */

    // Test basic equivalence:
    res &= rc::check(m_prefix + "Equivalence of View to model expression",
                     m_gen.generate(comptests::test_equivalence, 10. * eps));

    // Read-only element access
    res &= rc::check(m_prefix + "Element access via () and []",
                     m_gen.generate(comptests::test_element_access));
    res &= rc::check(m_prefix + "Element access via extract_block",
                     m_gen.generate(comptests::test_extract_block));
    res &= rc::check(m_prefix + "Data access via add_block_to",
                     m_gen.generate(comptests::test_add_block_to));
    res &= rc::check(m_prefix + "Read-only iterator of small matrices",
                     m_gen.generate(comptests::test_readonly_iterator));

    /* TODO enable if read-write views are implemented
    if (View::view_of_stored_matrix) {
        // Read-write element access
        res &= rc::check(
              m_prefix + "Altering elements via ()",
              m_gen.generate(comptests::test_setting_elements_indexed));
        res &= rc::check(
              m_prefix + "Altering elements via []",
              m_gen.generate(comptests::test_setting_elements_vectorised));
        res &= rc::check(m_prefix + "Altering elements via iterator",
                         m_gen.generate(comptests::test_readwrite_iterator));
    }
    */

    // Operations
    typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
          lazy_matrix_type;

    res &= rc::check(m_prefix + "Multiplication by scalar",
                     m_gen.generate(comptests::test_mutiply_scalar));
    res &= rc::check(m_prefix + "Divide by scalar",
                     m_gen.generate(comptests::test_divide_scalar));
    res &= rc::check(
          m_prefix + "Add a stored matrix",
          m_gen.generate(comptests::template test_add<stored_matrix_type>));
    res &= rc::check(
          m_prefix + "Add a lazy matrix",
          m_gen.generate(comptests::template test_add<lazy_matrix_type>));
    res &= rc::check(
          m_prefix + "Subtract a stored matrix",
          m_gen.generate(
                comptests::template test_subtract<stored_matrix_type>));
    res &= rc::check(
          m_prefix + "Multiply a stored matrix",
          m_gen.generate(
                comptests::template test_multiply_by<stored_matrix_type>));
    res &= rc::check(
          m_prefix + "Multiply a lazy matrix",
          m_gen.generate(
                comptests::template test_multiply_by<lazy_matrix_type>));

    return res;
}

}  // namespace view_tests
}  // namespace tests
}  // namespace linalgwrap
