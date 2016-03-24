#pragma once
#include "lazy_matrix_tests.hh"
#include "ScaleView.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

/* Namespace for default tests for views */
namespace view_tests {

/** Testing library for views.
 *
 * \tparam View        The view to test
 * \tparam ViewGenArg  The argument(s) needed to generate a view or
 *                     the equivalent model. Use a std::tuple or
 *                     similar if more than one argument is required
 * */
template <typename View, typename ViewGenArg>
class TestingLibrary
      : public lazy_matrix_tests::TestingLibrary<View, ViewGenArg> {
  public:
    typedef View view_type;
    typedef lazy_matrix_tests::TestingLibrary<View, ViewGenArg> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef ViewGenArg genarg_type;

    /** Construct a testing library instance in order to check
     *  all basic functionality of a lazy matrix.
     *
     *  \param args_generator   Function that generates the args for both
     *                          the following generators
     *  \param lazy_generator   Generator that takes the args and returns a lazy
     *                          matrix.
     *  \param model_generator  Generator that takes the args and returns a
     *                          model matrix
     *  \param prefix A prefix to use in the rapidcheck description string
     *  \param tolerance Numeric tolerance used for comparison in (some) checks.
     *                   It is not used for the check of basic equivalence,
     *                   where we require an agreement of 10*machine_epsilon.
     */
    TestingLibrary(
          std::function<genarg_type(void)> args_generator,
          std::function<view_type(genarg_type)> lazy_generator,
          std::function<stored_matrix_type(genarg_type)> model_generator,
          std::string prefix = "",
          double tolerance = TestConstants::default_num_tol);

    void run_checks() const;
};

/** Class to supply standard view_generator functionals for the TestingLibrary
 */
template <typename ViewTypes, typename MakeViewExtraArg>
struct StandardViewGenerators {
    // The types we use
    typedef typename ViewTypes::stored_matrix_type stored_matrix_type;

    // The derived default views:
    typedef typename ViewTypes::view_of_stored_type view_of_stored_type;
    typedef typename ViewTypes::view_of_scaleview_type view_of_scaleview_type;
    typedef typename ViewTypes::view_of_lazy_type view_of_lazy_type;

    // The argument and makeview function:
    typedef MakeViewExtraArg makeview_fctn_arg_type;

    /** Generate a view of a stored matrix */
    struct stored_view_generator {
        typedef typename view_of_stored_type::inner_matrix_type
              stored_matrix_type;
        typedef std::function<view_of_stored_type(
              stored_matrix_type&, makeview_fctn_arg_type)> makeview_fctn_type;

        view_of_stored_type operator()(
              std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args);

        stored_view_generator(makeview_fctn_type makeview_fctn);

      private:
        makeview_fctn_type m_makeview_fctn;
        std::shared_ptr<stored_matrix_type> m_stored_ptr;
    };

    /** Generate a view of a lazy matrix */
    struct lazy_view_generator {
        typedef typename view_of_lazy_type::inner_matrix_type lazy_matrix_type;
        typedef std::function<view_of_lazy_type(
              lazy_matrix_type&, makeview_fctn_arg_type)> makeview_fctn_type;

        view_of_lazy_type operator()(
              std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args);

        lazy_view_generator(makeview_fctn_type makeview_fctn);

      private:
        makeview_fctn_type m_makeview_fctn;
        std::shared_ptr<lazy_matrix_type> m_lazy_ptr;
    };

    /** Generate a view of a view of a stored matrix */
    struct view_view_generator {
        typedef typename view_of_scaleview_type::inner_matrix_type view_type;
        typedef std::function<view_of_scaleview_type(
              view_type&, makeview_fctn_arg_type)> makeview_fctn_type;

        static_assert(
              std::is_same<
                    view::ScaleView<typename view_type::inner_matrix_type>,
                    view_type>::value,
              "view_view_generator can only generate ScaleViews");

        view_of_scaleview_type operator()(
              std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args);

        view_view_generator(makeview_fctn_type makeview_fctn);

      private:
        makeview_fctn_type m_makeview_fctn;
        std::shared_ptr<stored_matrix_type> m_stored_ptr;
        std::shared_ptr<view_type> m_view_ptr;
    };
};

//
// --------------------------------------------------------
//

template <typename View, typename ViewGenArg>
TestingLibrary<View, ViewGenArg>::TestingLibrary(
      std::function<genarg_type(void)> args_generator,
      std::function<view_type(genarg_type)> view_generator,
      std::function<stored_matrix_type(genarg_type)> model_generator,
      std::string prefix, double tolerance)
      : base_type{args_generator, view_generator, model_generator, prefix,
                  tolerance} {}

template <typename View, typename ViewGenArg>
void TestingLibrary<View, ViewGenArg>::run_checks() const {

    // Run checks for the lazy part:
    base_type::run_checks();

    /* TODO enable if read-write views are implemented
    //Run the extras if we got a read-write view
    if (!View::is_const_view) {
        const double eps = std::numeric_limits<scalar_type>::epsilon();

        // Check copying
        CHECK(rc::check(m_prefix + "Test copying",
                         m_gen.generate(comptests::test_copy,eps)));

        // Read-write element access
        CHECK(rc::check(
              m_prefix + "Altering elements via ()",
              m_gen.generate(comptests::test_setting_elements_indexed)));
        CHECK(rc::check(
              m_prefix + "Altering elements via []",
              m_gen.generate(comptests::test_setting_elements_vectorised)));
        CHECK(rc::check(m_prefix + "Altering elements via iterator",
                         m_gen.generate(comptests::test_readwrite_iterator)));
    }
    */
}

//
// StandardViewGenerators
//
template <typename ViewTypes, typename MakeViewExtraArg>
typename StandardViewGenerators<ViewTypes,
                                MakeViewExtraArg>::view_of_stored_type
      StandardViewGenerators<ViewTypes,
                             MakeViewExtraArg>::stored_view_generator::
      operator()(
            std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args) {
    // Store the stored matrix in permanent memory:
    m_stored_ptr.reset(new stored_matrix_type{std::move(mat_args.first)});

    // Return a view to it:
    return m_makeview_fctn(std::ref(*m_stored_ptr), mat_args.second);
}

template <typename ViewTypes, typename MakeViewExtraArg>
StandardViewGenerators<ViewTypes, MakeViewExtraArg>::stored_view_generator::
      stored_view_generator(makeview_fctn_type makeview_fctn)
      : m_makeview_fctn(makeview_fctn) {}

template <typename ViewTypes, typename MakeViewExtraArg>
typename StandardViewGenerators<ViewTypes, MakeViewExtraArg>::view_of_lazy_type
      StandardViewGenerators<ViewTypes, MakeViewExtraArg>::lazy_view_generator::
      operator()(
            std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args) {
    // Store the lazy matrix in permanent memory:
    m_lazy_ptr.reset(new lazy_matrix_type{std::move(mat_args.first)});

    // Return a view to it:
    return m_makeview_fctn(std::ref(*m_lazy_ptr), mat_args.second);
}

template <typename ViewTypes, typename MakeViewExtraArg>
StandardViewGenerators<ViewTypes, MakeViewExtraArg>::lazy_view_generator::
      lazy_view_generator(makeview_fctn_type makeview_fctn)
      : m_makeview_fctn(makeview_fctn) {}

template <typename ViewTypes, typename MakeViewExtraArg>
typename StandardViewGenerators<ViewTypes,
                                MakeViewExtraArg>::view_of_scaleview_type
      StandardViewGenerators<ViewTypes, MakeViewExtraArg>::view_view_generator::
      operator()(
            std::pair<stored_matrix_type, makeview_fctn_arg_type> mat_args) {
    m_view_ptr.reset();

    // Store the stored matrix in permanent memory:
    m_stored_ptr.reset(new stored_matrix_type{std::move(mat_args.first)});

    // Construct a scale view and store it into the pointer:
    m_view_ptr.reset(new view_type{*m_stored_ptr, 1.0});

    // Return a view to the latter:
    return m_makeview_fctn(std::ref(*m_view_ptr), mat_args.second);
}

template <typename ViewTypes, typename MakeViewExtraArg>
StandardViewGenerators<ViewTypes, MakeViewExtraArg>::view_view_generator::
      view_view_generator(makeview_fctn_type makeview_fctn)
      : m_makeview_fctn(makeview_fctn) {}

}  // namespace view_tests
}  // namespace tests
}  // namespace linalgwrap
