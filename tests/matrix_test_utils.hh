#pragma once
#include <rapidcheck.h>
#include <rapidcheck/state.h>
#include <SmallMatrix.hh>
#include <list>

// Generators for extra matrices
#include "generators.hh"

// Numerical equality and comparison
#include "NumComp.hh"

// have an extra verbose output for rapidcheck function tests:
//#define HAVE_MATRIX_RC_CLASSIFY

// have a debug print of all generated matrices
// #define HAVE_MATRIX_DEBUG_PRINT

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for the components for the
 * Matrix model and standard matrix tests
 **/
namespace matrix_test_utils {

/** form the matrix product out = lhs*rhs
 */
template <typename Matrix1, typename Matrix2, typename Matrix3>
void matrix_product(const Matrix1& lhs, const Matrix2& rhs, Matrix3& out) {
    typedef typename Matrix1::size_type size_type;
    typedef typename Matrix1::scalar_type scalar_type;

    static_assert(std::is_same<typename Matrix1::scalar_type,
                               typename Matrix2::scalar_type>::value,
                  "All matrices need to have the same scalar type.");
    static_assert(std::is_same<typename Matrix2::scalar_type,
                               typename Matrix3::scalar_type>::value,
                  "All matrices need to have the same scalar type.");

    static_assert(std::is_same<typename Matrix1::size_type,
                               typename Matrix2::size_type>::value,
                  "All matrices need to have the same size type.");
    static_assert(std::is_same<typename Matrix2::size_type,
                               typename Matrix3::size_type>::value,
                  "All matrices need to have the same size type.");

    assert_size(lhs.n_rows(), out.n_rows());
    assert_size(rhs.n_cols(), out.n_cols());
    assert_size(lhs.n_cols(), rhs.n_rows());

    // TODO use iterator here?
    for (size_type i = 0; i < lhs.n_rows(); ++i) {
        for (size_type j = 0; j < rhs.n_cols(); ++j) {
            scalar_type sum = Constants<scalar_type>::zero;
            for (size_type k = 0; k < rhs.n_rows(); ++k) {
                sum += lhs(i, k) * rhs(k, j);
            }
            out(i, j) = sum;
        }
    }
}

template <typename Matrix>
struct MatrixSystemUnderTest {
    typedef Matrix matrix_type;

    //! The matrix we want to test
    matrix_type matrix;

    /** \brief construct a Test System to test a Matrix */
    MatrixSystemUnderTest(matrix_type mat) : matrix{mat} {}

    // Copy no matrices to internal storage.
    template <typename OtherMatrix>
    OtherMatrix& copy_to_internal_storage(OtherMatrix& m) {
        return m;
    }
};

template <typename MatrixModel, typename MatrixUnderTest>
struct TestingPolicy {
    typedef MatrixModel model_type;
    typedef MatrixSystemUnderTest<MatrixUnderTest> sut_type;

    /** Run some extra tests in order to compare/test model and sut.
     *
     * e.g. Test whethere a function in sut returns the expected result
     *versus
     * the model.
     *
     * The model and the sut are in the same state, i.e. should be
     *equivalent.
     * */
    static void extra_tests(const model_type&, const sut_type&) {
        // Do no extra tests
    }
};

template <typename TestingPolicy>
struct CommandBase : rc::state::Command<typename TestingPolicy::model_type,
                                        typename TestingPolicy::sut_type> {
  protected:
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;

    /** Perform the common checks for each command
     *
     * \param   model The model in the the state *before* the command has been
     *          applied i.e. the same state as run gets it.
     * \param   sut   The system under test in the state *after* the command has
     *          been applied i.e. after all the stuff of a command has been done
     *          to sut already.
     * */
    void run_common_tests(const model_type& model, const sut_type& sut) const {
#ifdef HAVE_MATRIX_DEBUG_PRINT
        std::cout << model << std::endl;
        std::cout << "------------------------------------" << std::endl;
#endif

        // Apply to a copy of the model
        model_type model_copy{model};
        this->apply(model_copy);

        // Check that the results are equivalent
        RC_ASSERT(NumComp::is_equal_matrix(model_copy, sut.matrix));

        // Do the additional tests requested:
        testpol_type::extra_tests(model_copy, sut);
    }
};

/** Add a Matrix to a MatrixUnderTest and to a Model.
 * Compare that the same thing arises.
 * Uses only low-level element-by-element operations
 * in the model, so almost all Matrix types should be
 * usable for MatrixModel or MatrixTermType */
template <typename TestingPolicy, typename MatrixTermType>
struct AddMatrix : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;
    typedef typename sut_type::stored_matrix_type stored_matrix_type;

    MatrixTermType term;

    AddMatrix(const model_type& model)
          : term{*FixedSizeMatrix<MatrixTermType>::fixed_size(
                        model.n_rows(), model.n_cols()).as("Matrix to add")} {}

    void apply(model_type& model) const override {
        // Make sure the size fits:
        RC_PRE(term.n_rows() == model.n_rows());
        RC_PRE(term.n_cols() == model.n_cols());

        // Add the stuff:
        // TODO use iterator here
        for (size_type i = 0; i < model.n_rows(); ++i) {
            for (size_type j = 0; j < model.n_cols(); ++j) {
                model(i, j) += term(i, j);
            }
        }
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Copy the matrix to internal storage if neccessary and
        // apply to sut
        const MatrixTermType& term_ref = sut.copy_to_internal_storage(term);
        sut.matrix = sut.matrix + term_ref;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "AddMatrix");
#endif

        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override {
        os << "AddMatrix (" << term << ")";
    }
};  // AddMatrix

/** Subtract a Matrix from a MatrixUnderTest and from a Model.
 * Compare that the same thing arises.
 * Uses only low-level element-by-element operations
 * in the model, so almost all Matrix types should be
 * usable for MatrixModel or MatrixTermType */
template <typename TestingPolicy, typename MatrixTermType>
struct SubtractMatrix : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;

    MatrixTermType term;

    SubtractMatrix(const model_type& model)
          : term{*FixedSizeMatrix<MatrixTermType>::fixed_size(model.n_rows(),
                                                              model.n_cols())
                        .as("Matrix to subtract.")} {}

    void apply(model_type& model) const override {
        // Make sure the size fits:
        RC_PRE(term.n_rows() == model.n_rows());
        RC_PRE(term.n_cols() == model.n_cols());

        // Add the stuff:
        // TODO use iterator here
        for (size_type i = 0; i < model.n_rows(); ++i) {
            for (size_type j = 0; j < model.n_cols(); ++j) {
                model(i, j) -= term(i, j);
            }
        }
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Copy the matrix to internal storage if neccessary and
        // apply to sut
        const MatrixTermType& term_ref = sut.copy_to_internal_storage(term);
        sut.matrix = sut.matrix - term_ref;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "SubtractMatrix");
#endif
        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override {
        os << "SubtractMatrix (" << term << ")";
    }
};  // SubtractMatrix

/** Perform a unary minus on the matrix, ie swap the signs
 * of all matrix entries. */
template <typename TestingPolicy>
struct UnaryMinusMatrix : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;

    UnaryMinusMatrix(const model_type&) {}

    void apply(model_type& model) const override {
        // swap the signs
        // TODO use iterator here
        for (size_type i = 0; i < model.n_rows(); ++i) {
            for (size_type j = 0; j < model.n_cols(); ++j) {
                model(i, j) = -model(i, j);
            }
        }
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Swap the signs:
        sut.matrix = -sut.matrix;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "UnaryMinusMatrix");
#endif
        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override { os << "UnaryMinusMatrix"; }
};  // UnaryMinusMatrix

/** Multiply a Matrix with a MatrixUnderTest and to a Model.
 * Compare that the same thing arises.
 * Uses only low-level element-by-element operations
 * in the model, so almost all Matrix types should be
 * usable for MatrixModel or MatrixTermType */
template <typename TestingPolicy, typename MatrixTermType>
struct MultiplyMatrix : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;

    MatrixTermType term;

    MultiplyMatrix(const model_type& model)
          : term{*FixedSizeMatrix<MatrixTermType>::fixed_size(
                        model.n_cols(),
                        *gen::inRange<size_type>(
                               1, TestConstants::max_matrix_size + 1)
                               .as("No of columns of multiplied matrix"))
                        .as("Matrix to multiply the current state with.")} {}

    void apply(model_type& model) const override {
        // Make sure the size fits:
        RC_PRE(model.n_cols() == term.n_rows());

        model_type out(model.n_rows(), term.n_cols(), false);
        matrix_product(model, term, out);
        model = out;
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Copy the matrix to internal storage if neccessary and
        // apply to sut
        const MatrixTermType& term_ref = sut.copy_to_internal_storage(term);
        sut.matrix = sut.matrix * term_ref;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "MultiplyMatrix");
#endif
        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override {
        os << "MultiplyMatrix (" << term << ")";
    }
};  // MultiplyMatrix

/** Multiply a Matrix by a scalar.
 * Compare that the same thing arises.
 */
template <typename TestingPolicy>
struct MultiplyScalar : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;

    scalar_type scalar;

    MultiplyScalar(const model_type&)
          : scalar{*rc::gen::nonZero<scalar_type>().as(
                  "Scalar to multiply with")} {}

    void apply(model_type& model) const override {
        // Multiply all entries
        // TODO use iterator here
        for (size_type i = 0; i < model.n_rows(); ++i) {
            for (size_type j = 0; j < model.n_cols(); ++j) {
                model(i, j) = scalar * model(i, j);
            }
        }
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Apply to sut:
        sut.matrix = scalar * sut.matrix;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "MultiplyScalar");
#endif
        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override {
        os << "MultiplyScalar (" << scalar << ")";
    }
};  // MultiplyScalar

/** Divide a Matrix by a scalar.
 * Compare that the same thing arises.
 */
template <typename TestingPolicy>
struct DivideScalar : CommandBase<TestingPolicy> {
    typedef CommandBase<TestingPolicy> base_type;
    typedef TestingPolicy testpol_type;
    typedef typename testpol_type::model_type model_type;
    typedef typename testpol_type::sut_type sut_type;
    typedef typename model_type::size_type size_type;
    typedef typename model_type::scalar_type scalar_type;

    scalar_type scalar;

    DivideScalar(const model_type&)
          : scalar{*rc::gen::nonZero<scalar_type>().as("Scalar to divide by")} {
    }

    void apply(model_type& model) const override {
        // Multiply all entries
        // TODO use iterator here
        for (size_type i = 0; i < model.n_rows(); ++i) {
            for (size_type j = 0; j < model.n_cols(); ++j) {
                model(i, j) = model(i, j) / scalar;
            }
        }
    }

    void run(const model_type& model, sut_type& sut) const override {
        // Apply to sut:
        sut.matrix = sut.matrix / scalar;

#ifdef HAVE_MATRIX_RC_CLASSIFY
        RC_CLASSIFY(true, "MultiplyScalar");
#endif
        base_type::run_common_tests(model, sut);
    }

    void show(std::ostream& os) const override {
        os << "MultiplyScalar (" << scalar << ")";
    }
};  // DivideScalar

}  // matrix_test_utils

}  // namespace tests
}  // namescpace linalgwrap
