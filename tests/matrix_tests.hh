#pragma once
#include <rapidcheck.h>
#include "NumComp.hh"
#include "Range.hh"
#include "generators.hh"

// have a debug print of all generated matrices
// #define HAVE_MATRIX_DEBUG_PRINT

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for the components for standard matrix tests.
 **/
namespace matrix_tests {
/** Low-level and basic function to form the matrix product out = lhs*rhs */
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

/** \brief Standard test functions which test a certain
 *  functionality by executing it in the SutMatrix
 *  and in a ModelMatrix and comparing the results
 *  afterwards.
 *
 *  \tparam CompMatrix  Model matrix used for comparison
 *  \tparam SutMatrix   System under test matrix.
 *                      The thing we test.
 **/
template <typename CompMatrix, typename SutMatrix>
struct ComparativeTests {
    typedef SutMatrix sutmat_type;
    typedef CompMatrix compmat_type;
    typedef typename SutMatrix::size_type size_type;
    typedef typename SutMatrix::scalar_type scalar_type;

    static_assert(
          std::is_same<size_type, typename CompMatrix::size_type>::value,
          "The size types of SutMatrix and CompMatrix have to agree");

    // TODO test swap function!

    /** Test read-only element access via () and [] at
     *  random places. Compare resulting values of the
     *  sut with the model */
    static void test_element_access(const compmat_type& model,
                                    const sutmat_type& sut);

    /** Test whether extracting a random block in the sut yields the data
     *  of that same block in the model i
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against values of the matrix itself.
     *  */
    static void test_extract_block(const compmat_type& model,
                                   const sutmat_type& sut);

    /** Test adding a block of the sut matrix to a different stored
     *  matrix. The yielded result is compared against the equivalent
     *  operation performed with the model.
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against values of the operation performed
     *       differently on the matrix itself.
     *  */
    static void test_add_block_to(const compmat_type& model,
                                  const sutmat_type& sut);

    /** Test whether using an iterator yields the same values as
     *  the entries in the model
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against the entries of the matrix.
     *  */
    static void test_iterator(const compmat_type& model,
                              const sutmat_type& sut);

    /** Test whether multiplication by a scalar yields the same
     *  values in model and sut */
    static void test_mutiply_scalar(const compmat_type& model,
                                    const sutmat_type& sut);

    /** Test whether division by a scalar yields the same
     *  values in model and sut */
    static void test_divide_scalar(const compmat_type& model,
                                   const sutmat_type& sut);

    /** Test whether addition of another arbitrary matrix gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherMatrix>
    static void test_add(const compmat_type& model, const sutmat_type& sut);

    /** Test whether subtraction of another arbitrary matrix gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherMatrix>
    static void test_subtract(const compmat_type& model,
                              const sutmat_type& sut);

    /** Test whether the multiplication with another arbitrary matrix
     *  gives rise to the same result in model and sut.
     */
    template <typename OtherMatrix>
    static void test_multiply_by(const compmat_type& model,
                                 const sutmat_type& sut);
};

//
// -------------------------------------------------------
//
template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_element_access(
      const compmat_type& model, const sutmat_type& sut) {
    size_type row = *gen::inRange<size_type>(0, model.n_rows()).as("row index");
    size_type col = *gen::inRange<size_type>(0, model.n_cols()).as("col index");
    size_type i = *gen::inRange<size_type>(0, model.n_rows() * model.n_cols())
                         .as("vectorised index");

    RC_ASSERT(NumComp::is_equal(model(row, col), sut(row, col)));
    RC_ASSERT(NumComp::is_equal(model[i], sut[i]));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_extract_block(
      const compmat_type& model, const sutmat_type& sut) {
    Range<size_type> rows =
          *gen::scale(2.0, gen::range_within<size_type>(0, model.n_rows()))
                 .as("extract_block row range");
    Range<size_type> cols =
          *gen::scale(2.0, gen::range_within<size_type>(0, model.n_cols()))
                 .as("extract_block col range");

    // Extract a block from the matrix:
    // The type is determined by the return value of the function.
    // It is pretty sure a matrix, so the next block of code works.
    auto block = sut.extract_block(rows, cols);

    // check that it agrees with the model:
    for (size_type row : rows) {
        for (size_type col : cols) {
            RC_ASSERT(NumComp::is_equal(model(row, col), sut(row, col)));
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_add_block_to(
      const compmat_type& model, const sutmat_type& sut) {
    // the return type of extract_block and the type add_block_to expects
    // is identical. Hence we can use this trick in order to obtain the
    // type of matrix we need to generate.
    typedef decltype(sut.extract_block(
          Range<size_type>{0, 0}, Range<size_type>{0, 0})) block_matrix_type;

    // Then generate two sizes:
    auto n_rows = *gen::inRange<size_type>(0, model.n_rows())
                         .as("Number of rows of the matrix to add a block to");
    auto n_cols = *gen::inRange<size_type>(0, model.n_cols()).as(
          "Number of columns of the matrix to add a block to");

    // Generate the offsets:
    auto r_offset = *gen::inRange<size_type>(0, model.n_rows() - n_rows)
                           .as("Start row for add_block_to");
    auto c_offset = *gen::inRange<size_type>(0, model.n_cols() - n_cols)
                           .as("Start col for add_block_to");

    // Generate an arbitrary factor:
    auto c_this = *gen::arbitrary<scalar_type>().as("c_this");

    // Generate an arbitrary matrix of this size:
    block_matrix_type inmat =
          *gen::fixed_size<block_matrix_type>(n_rows, n_cols);

    // Copy the original
    block_matrix_type inmat_copy{inmat};

    // Perform the add_block_to operation
    sut.add_block_to(inmat, r_offset, c_offset, c_this);

    // Check that the result is correct:
    for (auto row : range(n_rows)) {
        for (auto col : range(n_cols)) {
            RC_ASSERT(NumComp::is_equal(
                  inmat(row, col),
                  inmat_copy(row, col) +
                        c_this * sut(row + r_offset, col + c_offset)));
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_iterator(
      const compmat_type& model, const sutmat_type& sut) {
    auto it_const = sut.cbegin();
    auto it = sut.begin();

    for (size_type i = 0; i < model.n_rows(); ++i) {
        for (size_type j = 0; j < model.n_cols(); ++j) {
            // Assert that const and non-const iterators
            // behave identically:
            RC_ASSERT(*it == *it_const);
            RC_ASSERT(it.row() == it_const.row());
            RC_ASSERT(it.col() == it_const.col());

            if (!NumComp::is_equal(model(i, j), 0.,
                                   0.01 * TestConstants::default_num_tol,
                                   false)) {
                // Larger than zero element -> expect to be present
                RC_ASSERT(it.row() == i);
                RC_ASSERT(it.col() == j);
                RC_ASSERT(NumComp::is_equal(*it, model(i, j)));
            }

            if (it.indices() <= std::make_pair(i, j)) {
                // increment iterators if neccessary and continue
                ++it;
                ++it_const;
            }
        }
    }

    for (auto it = std::begin(sut); it != std::end(sut); ++it) {
        RC_ASSERT(NumComp::is_equal(*it, model(it.row(), it.col())));
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_mutiply_scalar(
      const compmat_type& model, const sutmat_type& sut) {
    // Generate an arbitrary factor:
    auto c = *gen::arbitrary<scalar_type>().as("Coefficient");

    // Do the multiplication:
    auto res = sut * c;
    auto res2 = c * sut;

    compmat_type res_model{model.n_rows(), model.n_cols(), false};
    for (auto row : range(model.n_rows())) {
        for (auto col : range(model.n_cols())) {
            res_model(row, col) = model(row, col) * c;
        }
    }

    RC_ASSERT(NumComp::is_equal_matrix(res, res_model));
    RC_ASSERT(NumComp::is_equal_matrix(res2, res_model));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_divide_scalar(
      const compmat_type& model, const sutmat_type& sut) {
    // Generate an arbitrary factor:
    auto c = *gen::nonZero<scalar_type>().as("Coefficient");

    // Do the multiplication:
    auto res = sut / c;

    compmat_type res_model{model.n_rows(), model.n_cols(), false};
    for (auto row : range(model.n_rows())) {
        for (auto col : range(model.n_cols())) {
            res_model(row, col) = model(row, col) / c;
        }
    }
    RC_ASSERT(NumComp::is_equal_matrix(res, res_model));
}

template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_multiply_by(
      const compmat_type& model, const sutmat_type& sut) {
    // The size of the other matrix to multiply mlhs with:
    // Note: the RHS of inRange is exclusive
    auto othersize =
          *gen::inRange<size_type>(1, TestConstants::max_matrix_size + 1)
                 .as("Number of columns of the rhs matrix");

    // Generate another matrix:
    OtherMatrix mrhs = *gen::fixed_size<OtherMatrix>(model.n_cols(), othersize)
                              .as("Multiplication rhs");

    // Do the low-level multiplication:
    compmat_type result_model(model.n_rows(), othersize, false);
    matrix_product(model, mrhs, result_model);

    // multipy them
    auto result = sut * mrhs;

    // Check that the results are equivalent
    RC_ASSERT(NumComp::is_equal_matrix(result_model, result));
}

/** Test whether addition of another arbitrary matrix gives
 *  rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_add(
      const compmat_type& model, const sutmat_type& sut) {
    // generate another matrix of the same size:
    auto madd = *gen::fixed_size<OtherMatrix>(model.n_rows(), model.n_cols())
                       .as("Matrix to add");

    // Perform the operation
    auto res = sut + madd;

    // and on the model:
    compmat_type res_model(model.n_rows(), model.n_cols(), false);
    for (auto row : range(model.n_rows())) {
        for (auto col : range(model.n_cols())) {
            res_model(row, col) = model(row, col) + madd(row, col);
        }
    }

    // Check that the results are equivalent
    RC_ASSERT(NumComp::is_equal_matrix(res_model, res));
}

/** Test whether subtraction of another arbitrary matrix gives
 *  rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_subtract(
      const compmat_type& model, const sutmat_type& sut) {
    // generate another matrix of the same size:
    auto msub = *gen::fixed_size<OtherMatrix>(model.n_rows(), model.n_cols())
                       .as("Matrix to subtract");

    // Perform the operation
    auto res = sut - msub;

    // and on the model:
    compmat_type res_model(model.n_rows(), model.n_cols(), false);
    for (auto row : range(model.n_rows())) {
        for (auto col : range(model.n_cols())) {
            res_model(row, col) = model(row, col) - msub(row, col);
        }
    }

    // Check that the results are equivalent
    RC_ASSERT(NumComp::is_equal_matrix(res_model, res));
}

}  // matrix_tests
}  // namespace tests
}  // namescpace linalgwrap
