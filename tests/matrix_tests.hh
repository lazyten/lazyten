#pragma once
#include <rapidcheck.h>
#include "NumComp.hh"
#include "Range.hh"
#include "generators.hh"
#include "TestConstants.hh"
#include <functional>

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
void matrix_product(const Matrix1& lhs, const Matrix2& rhs, Matrix3& out);

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

    /** std::function type which describes the call signature of all
     *  static test functions in this class. */
    typedef std::function<void(const compmat_type&, const sutmat_type&,
                               const double)> testfunction_type;

    static_assert(
          std::is_same<size_type, typename CompMatrix::size_type>::value,
          "The size types of SutMatrix and CompMatrix have to agree");

    // TODO test swap function!

    /** Test wheather the two matrices are identical */
    static void test_equivalence(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test copying the sut */
    static void test_copy(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test read-only element access via () and [] at
     *  random places. Compare resulting values of the
     *  sut with the model */
    static void test_element_access(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test read-write element access via () at a random place.
     *  compare result against a model at all places except the
     *  changed entry.
     */
    static void test_setting_elements_indexed(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test read-write element access via [] at a random place.
     *  compare result against a model at all places except the
     *  changed entry.
     */
    static void test_setting_elements_vectorised(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether extracting a random block in the sut yields the data
     *  of that same block in the model i
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against values of the matrix itself.
     *  */
    static void test_extract_block(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test adding a block of the sut matrix to a different stored
     *  matrix. The yielded result is compared against the equivalent
     *  operation performed with the model.
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against values of the operation performed
     *       differently on the matrix itself.
     *  */
    static void test_add_block_to(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether using an iterator yields the same values as
     *  the entries in the model
     *
     * \note if sut and model are the same object this function compares
     *       the values yielded against the entries of the matrix.
     *  */
    static void test_readonly_iterator(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether modifying element entries using an iterator
     *  works.
     **/
    static void test_readwrite_iterator(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether multiplication by a scalar yields the same
     *  values in model and sut */
    static void test_mutiply_scalar(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether division by a scalar yields the same
     *  values in model and sut */
    static void test_divide_scalar(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether addition of another arbitrary matrix gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherMatrix>
    static void test_add(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether subtraction of another arbitrary matrix gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherMatrix>
    static void test_subtract(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** Test whether the multiplication with another arbitrary matrix
     *  gives rise to the same result in model and sut.
     */
    template <typename OtherMatrix>
    static void test_multiply_by(
          const compmat_type& model, const sutmat_type& sut,
          const double tolerance = TestConstants::default_num_tol);

    /** \brief Helper class to simplify setup and calling the above functions
     *  from rapidcheck.
     *
     * The idea is that we construct an object of this type and supply it with
     * the information needed to setup the model and the sut in order to
     * call the test functions above.
     *
     * On call of the generate function a std::function<void(void)> is
     * returned which calls the actual testfunction with a model and sut,
     * generated by the sut_generator and model_generator using common
     * arguments which are in turn generated by the args_generator.
     * See stored_matrix_tests.hh for an example.
     *
     * \tparam GenArg The argument type(s) we need to generate a sut matrix
     * or a model matrix (e.g. a matrix and a factor or a matrix
     * and some sizes ...). Use a std::tuple if more than one arg is required.
     */
    template <typename GenArg>
    class RapidcheckTestableGenerator {
      public:
        /** \brief Construct a RapidcheckCallGenerator object.
         *
         * \param sutgenerator  A function that derives a System-under-test
         *                      matrix from an randomly generated core model
         *                      matrix of type compmat_type.
         * \param modelgenerator  A function that derives a compmat_type model
         *                        matrix in the same state as the test matrix
         *                        yielded by sutgenerator when applied to the
         *                        same core model matrix.
         * \param tolerance     Default numeric tolerance to use.
         */
        RapidcheckTestableGenerator(
              std::function<GenArg(void)> args_generator,
              std::function<sutmat_type(GenArg)> sut_generator,
              std::function<compmat_type(GenArg)> model_generator,
              double tolerance = TestConstants::default_num_tol);

        /* Return a std::function object that calls the supplied
         * function func with the generated model and sut.
         *
         * Uses the tolerance supplied upon class construction
         */
        std::function<void(void)> generate(testfunction_type func) const;

        /* Return a std::function object that calls the supplied
         * function func with the generated model and sut.
         *
         * Uses a different tolerance.
         */
        std::function<void(void)> generate(testfunction_type func,
                                           double tolerance) const;

      private:
        std::function<GenArg(void)> m_args_generator;
        std::function<sutmat_type(GenArg)> m_sut_generator;
        std::function<compmat_type(GenArg)> m_model_generator;
        double m_tolerance;
    };
};

//
// -------------------------------------------------------
//

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

//
// ----------
//
//
template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_equivalence(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    RC_ASSERT(NumComp::is_equal_matrix(model, sut, tolerance));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_copy(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    sutmat_type copy{sut};

    // check that it is identical to the model
    RC_ASSERT(NumComp::is_equal_matrix(copy, model, tolerance));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_element_access(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    size_type row = *gen::inRange<size_type>(0, model.n_rows()).as("row index");
    size_type col = *gen::inRange<size_type>(0, model.n_cols()).as("col index");
    size_type i = *gen::inRange<size_type>(0, model.n_rows() * model.n_cols())
                         .as("vectorised index");

    RC_ASSERT(NumComp::is_equal(model(row, col), sut(row, col), tolerance));
    RC_ASSERT(NumComp::is_equal(model[i], sut[i], tolerance));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_setting_elements_indexed(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {

    auto modify_row =
          *gen::inRange<size_type>(0, model.n_rows()).as("Row to modify");
    auto modify_col =
          *gen::inRange<size_type>(0, model.n_cols()).as("Col to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    sutmat_type sut_copy{sut};

    // Modify the value
    sut_copy(modify_row, modify_col) = value;

    // Check it:
    for (auto row : range(model.n_rows())) {
        for (auto col : range(model.n_cols())) {
            if (row == modify_row && col == modify_col) {
                RC_ASSERT(
                      NumComp::is_equal(sut_copy(row, col), value, tolerance));
            } else {
                RC_ASSERT(NumComp::is_equal(sut_copy(row, col), model(row, col),
                                            tolerance));
            }
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_setting_elements_vectorised(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    auto index = *gen::inRange<size_type>(0, model.n_rows() * model.n_cols())
                        .as("Index to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    sutmat_type sut_copy{sut};

    // Modify the value
    sut_copy[index] = value;

    // Check it:
    for (auto i : range(model.n_rows() * model.n_cols())) {
        if (i == index) {
            RC_ASSERT(NumComp::is_equal(sut_copy[i], value, tolerance));
        } else {
            RC_ASSERT(NumComp::is_equal(sut_copy[i], model[i], tolerance));
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_extract_block(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {

    size_type min_range_size = 1;
    Range<size_type> rows =
          *gen::scale(2.0, gen::range_within<size_type>(0, model.n_rows(),
                                                        min_range_size))
                 .as("extract_block row range");
    Range<size_type> cols =
          *gen::scale(2.0, gen::range_within<size_type>(0, model.n_cols(),
                                                        min_range_size))
                 .as("extract_block col range");

    // Extract a block from the matrix:
    // The type is determined by the return value of the function.
    // It is pretty sure a matrix, so the next block of code works.
    auto block = sut.extract_block(rows, cols);

    // check that it agrees with the model:
    for (size_type row : rows) {
        for (size_type col : cols) {
            size_type i = row - rows.first();
            size_type j = col - cols.first();
            RC_ASSERT(
                  NumComp::is_equal(model(row, col), block(i, j), tolerance));
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_add_block_to(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    // the return type of extract_block and the type add_block_to expects
    // is identical. Hence we can use this trick in order to obtain the
    // type of matrix we need to generate.
    typedef decltype(sut.extract_block(
          Range<size_type>{0, 0}, Range<size_type>{0, 0})) block_matrix_type;

    // Then generate two sizes:
    size_type min_size = 1;
    auto n_rows = *gen::inRange<size_type>(min_size, model.n_rows() + 1)
                         .as("Number of rows of the matrix to add a block to");
    auto n_cols = *gen::inRange<size_type>(min_size, model.n_cols() + 1).as(
          "Number of columns of the matrix to add a block to");

    // Generate the offsets:
    auto r_offset = *gen::inRange<size_type>(0, model.n_rows() - n_rows + 1)
                           .as("Start row for add_block_to");
    auto c_offset = *gen::inRange<size_type>(0, model.n_cols() - n_cols + 1)
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
                        c_this * sut(row + r_offset, col + c_offset),
                  tolerance));
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_readonly_iterator(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    auto it_const = sut.cbegin();
    auto it = sut.begin();

    for (size_type i = 0; i < model.n_rows(); ++i) {
        for (size_type j = 0; j < model.n_cols(); ++j) {
            // Assert that const and non-const iterators
            // behave identically:
            RC_ASSERT(*it == *it_const);
            RC_ASSERT(it.row() == it_const.row());
            RC_ASSERT(it.col() == it_const.col());

            if (!NumComp::is_equal(
                      model(i, j), 0.,
                      10. * std::numeric_limits<scalar_type>::epsilon(),
                      false)) {
                // Larger than zero element -> expect to be present
                RC_ASSERT(it.row() == i);
                RC_ASSERT(it.col() == j);
                RC_ASSERT(NumComp::is_equal(*it, model(i, j), tolerance));
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
void ComparativeTests<CompMatrix, SutMatrix>::test_readwrite_iterator(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
    auto modify_row =
          *gen::inRange<size_type>(0, model.n_rows()).as("Row to modify");
    auto modify_col =
          *gen::inRange<size_type>(0, model.n_cols()).as("Col to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Modify an entry
    sutmat_type sut_copy{sut};
    for (auto it = sut_copy.begin(); it != sut_copy.end(); ++it) {
        if (it.row() == modify_row && it.col() == modify_col) {
            *it = value;
            break;
        }
    }

    for (size_type i = 0; i < model.n_rows(); ++i) {
        for (size_type j = 0; j < model.n_cols(); ++j) {
            if (i == modify_row && j == modify_col) {
                RC_ASSERT(NumComp::is_equal(sut_copy(i, j), value, tolerance));
            } else {
                RC_ASSERT(NumComp::is_equal(model(i, j), sut_copy(i, j),
                                            tolerance));
            }
        }
    }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_mutiply_scalar(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
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

    RC_ASSERT(NumComp::is_equal_matrix(res, res_model, tolerance));
    RC_ASSERT(NumComp::is_equal_matrix(res2, res_model, tolerance));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_divide_scalar(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
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
    RC_ASSERT(NumComp::is_equal_matrix(res, res_model, tolerance));
}

template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_multiply_by(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
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
    RC_ASSERT(NumComp::is_equal_matrix(result_model, result, tolerance));
}

/** Test whether addition of another arbitrary matrix gives
 *  rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_add(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
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
    RC_ASSERT(NumComp::is_equal_matrix(res_model, res, tolerance));
}

/** Test whether subtraction of another arbitrary matrix gives
 *  rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_subtract(
      const compmat_type& model, const sutmat_type& sut,
      const double tolerance) {
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
    RC_ASSERT(NumComp::is_equal_matrix(res_model, res, tolerance));
}

template <typename CompMatrix, typename SutMatrix>
template <typename GenArg>
ComparativeTests<CompMatrix, SutMatrix>::RapidcheckTestableGenerator<GenArg>::
      RapidcheckTestableGenerator(
            std::function<GenArg(void)> args_generator,
            std::function<sutmat_type(GenArg)> sut_generator,
            std::function<compmat_type(GenArg)> model_generator,
            double tolerance)
      : m_args_generator(args_generator),
        m_sut_generator(sut_generator),
        m_model_generator(model_generator),
        m_tolerance(tolerance) {}

template <typename CompMatrix, typename SutMatrix>
template <typename GenArg>
std::function<void(void)>
ComparativeTests<CompMatrix, SutMatrix>::RapidcheckTestableGenerator<
      GenArg>::generate(testfunction_type func, double tolerance) const {
    // Return a lamda which calls the test appropriately.
    return [=] {
        // Generate the args:
        GenArg arg = m_args_generator();

        // Generate the sut:
        sutmat_type s = m_sut_generator(arg);

        // Generate the equivalent model:
        compmat_type m = m_model_generator(arg);

        // Finally call the test.
        func(m, s, tolerance);
    };
}

template <typename CompMatrix, typename SutMatrix>
template <typename GenArg>
std::function<void(void)>
ComparativeTests<CompMatrix, SutMatrix>::RapidcheckTestableGenerator<
      GenArg>::generate(testfunction_type func) const {
    return generate(func, m_tolerance);
}

}  // matrix_tests
}  // namespace tests
}  // namescpace linalgwrap
