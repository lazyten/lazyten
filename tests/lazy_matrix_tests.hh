#pragma once
#include <rapidcheck.h>
#include <LazyMatrixExpression.hh>
#include "matrix_tests.hh"
#include "rapidcheck_utils.hh"

// have an extra verbose output for rapidcheck function tests:
//# define HAVE_LAZYMATRIX_RC_CLASSIFY

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for the components for the
 * LazyMatrix model and standard lazy matrix tests
 **/
namespace lazy_matrix_tests {

/** Default testing routines for lazy matrices to test
 *  the required functionality of a lazy matrix in the
 *  current state
 *
 *  Achieve this by comparing the outcome of the lazy
 *  matrix operation with the equivalent operation
 *  performed on a model
 *
 *  \tparam CompMatrix  Model matrix used for comparison
 *  \tparam SutMatrix   System under test matrix.
 *                      The thing we test.
 *  */
template <typename CompMatrix, typename SutMatrix>
struct FunctionalityTests
      : public matrix_tests::ComparativeTests<CompMatrix, SutMatrix> {
    typedef matrix_tests::ComparativeTests<CompMatrix, SutMatrix> base_type;
    typedef typename base_type::sutmat_type sutmat_type;
    typedef typename base_type::compmat_type compmat_type;
    typedef typename base_type::size_type size_type;
    typedef typename sutmat_type::stored_matrix_type stored_matrix_type;

    /** \brief Run all standard tests in order to compare/test model and sut. */
    static void run_all_tests(const compmat_type& model,
                              const sutmat_type& sut) {
        base_type::test_copy(model, sut);
        base_type::test_element_access(model, sut);
        base_type::test_equivalence(model, sut);
        base_type::test_extract_block(model, sut);
        base_type::test_add_block_to(model, sut);
        base_type::test_readonly_iterator(model, sut);
        base_type::template test_multiply_by<stored_matrix_type>(model, sut);
        test_convert_to_stored(model, sut);
    }

    /** Test conversion of the sut metrix to its stored_matrix_type */
    static void test_convert_to_stored(const compmat_type& model,
                                       const sutmat_type& sut) {
        stored_matrix_type sm = static_cast<stored_matrix_type>(sut);

        // Check that it is equivalent to the model:
        RC_ASSERT(NumComp::is_equal_matrix(sm, model));
    }
};

}  // namespace lazy_matrix_tests
}  // namespace tests
}  // namescpace linalgwrap
