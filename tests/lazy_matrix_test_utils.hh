#pragma once
#include <rapidcheck.h>
#include <rapidcheck/state.h>
#include <LazyMatrixExpression.hh>
#include "generators.hh"
#include "NumComp.hh"
#include "matrix_test_utils.hh"
#include "rapidcheck_utils.hh"

// have an extra verbose output for rapidcheck function tests:
//# define HAVE_LAZYMATRIX_RC_CLASSIFY

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for the components for the
 * LazyMatrix model and standard lazy matrix tests
 **/
namespace lazy_matrix_test_utils {

template <typename LazyMatrix>
struct LazyMatrixSystemUnderTest {
    typedef LazyMatrix matrix_type;
    typedef typename matrix_type::stored_matrix_type stored_matrix_type;

    //! The stored matrices, which are indirectly referred to by the
    //  matrix object below (due to subscriptions)
    //  We use a reference from global scope. This assures, that
    //  the objects are deleted after the actual matrices.
    std::list<stored_matrix_type>& stored_matrices;

    //! The matrix we want to test
    matrix_type matrix;

    /** \brief construct a Test System to test a Matrix
     *
     * \param mat                     The initial matrix
     * \param stored_matrices_storage An empty list of stored matrices which
     *              will be used in order to store the stored matrices we refer
     *              to in global scope.
     * */
    LazyMatrixSystemUnderTest(
          matrix_type mat,
          std::list<stored_matrix_type>& stored_matrices_storage)
          : stored_matrices{stored_matrices_storage}, matrix{mat} {}

    /* Offer the class a chance to copy a matrix into internal storage
     * We will not do so unless it is a stored matrix*/
    template <typename OtherMatrix>
    OtherMatrix& copy_to_internal_storage(OtherMatrix& m) {
        return m;
    }

    stored_matrix_type& copy_to_internal_storage(stored_matrix_type& m) {
        stored_matrices.push_back(m);
        return stored_matrices.back();
    }

    const stored_matrix_type& copy_to_internal_storage(
          const stored_matrix_type& m) {
        stored_matrices.push_back(m);
        return stored_matrices.back();
    }
};

template <typename MatrixModel, typename MatrixUnderTest>
struct LazyMatrixTestingPolicy {
    typedef MatrixModel model_type;
    typedef LazyMatrixSystemUnderTest<MatrixUnderTest> sut_type;
    typedef typename sut_type::matrix_type matrix_type;
    typedef typename matrix_type::stored_matrix_type stored_matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    /** Run some extra tests in order to compare/test model and sut.
     *
     * In this case we test the operator* with a Stored Matrix and
     * the fill() function.
     * */
    static void extra_tests(const model_type& model, const sut_type& sut) {
        test_mult_by_stored(model, sut);
        test_fill(model, sut);
    }

    static void test_mult_by_stored(const model_type& model,
                                    const sut_type& sut) {
        // The size of the other matrix to multiply with:
        // The RHS of inRange is exclusive
        auto othersize =
              *gen::inRange<size_type>(1, TestConstants::max_matrix_size + 1)
                     .as("Columns of the StoredMatrix to multiply");

        // Generate a random stored matrix
        auto mat = *FixedSizeMatrix<stored_matrix_type>::fixed_size(
                          model.n_cols(), othersize)
                          .as("Stored matrix to multiply with");

        stored_matrix_type result_model(model.n_rows(), othersize, false);
        matrix_test_utils::matrix_product(model, mat, result_model);

        stored_matrix_type result_matrix = sut.matrix * mat;

        // Check that the results are equivalent
        RC_ASSERT(NumComp::is_equal_matrix(result_model, result_matrix));
    }

    static void test_fill(const model_type& model, const sut_type& sut) {
        size_type start_row = *gen::inRange<size_type>(0, model.n_rows());
        size_type start_col = *gen::inRange<size_type>(0, model.n_cols());

        size_type n_rows =
              *gen::inRange<size_type>(1, model.n_rows() - start_row + 1);
        size_type n_cols =
              *gen::inRange<size_type>(1, model.n_cols() - start_col + 1);

        // Extract a block from the matrix:
        SmallMatrix<scalar_type> block(n_rows, n_cols, false);
        sut.matrix.fill(start_row, start_col, block);

        // check that it agrees with the model:
        for (size_type i = start_row; i < start_row + n_rows; ++i) {
            for (size_type j = start_col; j < start_col + n_cols; ++j) {
                RC_ASSERT(NumComp::is_equal(model(i, j), sut.matrix(i, j)));
            }
        }
    }
};

template <typename LazyMatrix>
struct CommandDefaultTypes {
    typedef LazyMatrix matrix_type;

    //! The default stored matrix type used in the operators
    typedef typename matrix_type::stored_matrix_type stored_matrix_type;

    //! The default lazy matrix type used in the operators
    typedef LazyMatrixWrapper<stored_matrix_type, stored_matrix_type>
          lazy_matrix_type;
};

template <typename ModelMatrix, typename LazyMatrix,
          typename CommandDefaultTypes = CommandDefaultTypes<LazyMatrix>>
struct TestLibrary {
    typedef LazyMatrix matrix_type;
    typedef ModelMatrix model_matrix_type;
    typedef typename matrix_type::stored_matrix_type stored_matrix_type;
    typedef CommandDefaultTypes command_default_types;

    //! The testing policy we will use for checking
    typedef lazy_matrix_test_utils::LazyMatrixTestingPolicy<
          model_matrix_type, matrix_type> testing_policy;

    //! The default commands defined
    typedef matrix_test_utils::MultiplyMatrix<
          testing_policy, typename command_default_types::lazy_matrix_type>
          op_MultiplyLazy;
    typedef matrix_test_utils::UnaryMinusMatrix<testing_policy> op_UnaryMinus;
    typedef matrix_test_utils::MultiplyScalar<testing_policy> op_MultScalar;
    typedef matrix_test_utils::DivideScalar<testing_policy> op_DivideScalar;
    typedef matrix_test_utils::AddMatrix<
          testing_policy, typename command_default_types::lazy_matrix_type>
          op_AddLazy;
    typedef matrix_test_utils::AddMatrix<
          testing_policy, typename command_default_types::stored_matrix_type>
          op_AddStored;
    typedef matrix_test_utils::SubtractMatrix<
          testing_policy, typename command_default_types::lazy_matrix_type>
          op_SubtractLazy;
    typedef matrix_test_utils::SubtractMatrix<
          testing_policy, typename command_default_types::stored_matrix_type>
          op_SubtractStored;

    template <typename GenFunc>
    void run_check(const model_matrix_type& initial_state,
                   GenFunc&& generation_func, double scale = 1.0) const {

        //! The global space used for stroing matrices the lazy matrix
        // expressions
        //  point to indirectly
        std::list<stored_matrix_type> stored_matrices;

        // The resulting System Under Test type:
        typedef typename testing_policy::sut_type sut_type;

        // Setup the initial system and the initial_sut:
        LazyMatrixWrapper<stored_matrix_type, model_matrix_type> wrap{
              initial_state};
        matrix_type model{wrap};
        sut_type initial_sut{model, stored_matrices};

        // TODO: Until something better exists
        // Run it through rapidcheck
        rapidcheck_utils::state_check_scaled(initial_state, initial_sut, scale,
                                             generation_func);
    }
};

}  // namespace lazy_matrix_tests
}  // namespace tests
}  // namescpace linalgwrap
