//
// Copyright (C) 2016-17 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#pragma once
#include "RCTestableGenerator.hh"
#include "generators.hh"
#include "indexable_tests.hh"
#include "macro_defs.hh"
#include "rapidcheck_utils.hh"
#include "vector_tests.hh"
#include <functional>
#include <lazyten/TestingUtils.hh>
#include <lazyten/TypeUtils.hh>

namespace lazyten {
namespace tests {
using namespace rc;
using namespace krims;

/** Namespace for the components for standard matrix tests.
 **/
namespace matrix_tests {
/** Low-level and basic function to form the matrix product out = lhs*rhs */
template <typename Matrix1, typename Matrix2, typename Matrix3>
void matrix_product(const Matrix1& lhs, const Matrix2& rhs, Matrix3& out);

/** Use a really stupid algorithm to return the transpose of the given matrix
 * inside out)
 */
template <typename Matrix1, typename Matrix2>
void matrix_transpose(const Matrix1& in, Matrix2& trans);

/** Low-level and basic function to apply a matrix to a vector */
template <typename Matrix, typename Vector1, typename Vector2>
void matrix_apply(const Matrix& lhs, const Vector1& rhs, Vector2& out);

/** \brief Standard test functions which test a certain
 *  functionality by executing it in the Sut indexable
 *  and in a Model matrix and comparing the results
 *  afterwards.
 *
 *  \tparam Model  Model matrix used for comparison
 *  \tparam Sut   System under test matrix.
 *                      The thing we test.
 **/
template <typename Model, typename Sut>
struct ComparativeTests {
  typedef Model model_type;
  typedef Sut sut_type;
  typedef typename sut_type::size_type size_type;
  typedef typename sut_type::scalar_type scalar_type;

  static_assert(std::is_same<size_type, typename Model::size_type>::value,
                "The size types of Sut and Model have to agree");

  // The related comptest types for indexables and vectors
  typedef indexable_tests::ComparativeTests<Model, Sut> idx_cmptests;

  // TODO test swap function!

  // TODO use stuff from indexable_tests.hh and vector_tests.hh
  // TODO missing from indexable_tests
  //  - test_dot

  // TODO temporary
  typedef model_type compmat_type;
  typedef sut_type sutmat_type;

  /** Test whether the two matrices are identical */
  static void test_equivalence(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test copying the sut */
  lazyten_declare_comptest(test_copy);

  /** Test read-only element access via () and [] at
   *  random places. Compare resulting values of the
   *  sut with the model */
  static void test_element_access(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test read-write element access via () at a random place.
   *  compare result against a model at all places except the
   *  changed entry.
   */
  static void test_setting_elements_indexed(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test read-write element access via [] at a random place.
   *  compare result against a model at all places except the
   *  changed entry.
   */
  static void test_setting_elements_vectorised(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether extracting a random block in the sut yields the data
   *  of that same block in the model i
   *
   * \note if sut and model are the same object this function compares
   *       the values yielded against values of the matrix itself.
   *  */
  lazyten_declare_comptest(test_extract_block);

  /** Test whether extracting a random block in the sut yields the data
   *  of that same block in the model i
   *  We use transpose operation mode here
   *
   * \note if sut and model are the same object this function compares
   *       the values yielded against values of the matrix itself.
   *  */
  lazyten_declare_comptest(test_extract_transpose_block);

  // TODO override version of indexable_tests.hh
  /** Test whether using an iterator yields the same values as
   *  the entries in the model
   *
   * \note if sut and model are the same object this function compares
   *       the values yielded against the entries of the matrix.
   *  */
  static void test_readonly_iterator(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether modifying element entries using an iterator
   *  works.
   **/
  static void test_readwrite_iterator(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test the l1 norm function */
  lazyten_declare_comptest(test_norm_l1);

  /** Test the linf norm function */
  lazyten_declare_comptest(test_norm_linf);

  /** Test the frobenius norm functions (norm_frobenius and
   * norm_frobenius_squared) */
  lazyten_declare_comptest(test_norm_frobenius);

  /** Test the trace function */
  lazyten_declare_comptest(test_trace);

  /** Test the functions min and max*/
  lazyten_declare_comptest(test_minmax);

  /** Test the elementwise functions abs, conj, sqrt and square */
  lazyten_declare_comptest(test_elementwise);

  /** Test the accumulate function */
  lazyten_declare_comptest(test_accumulate);

  /** Test whether multiplication by a scalar yields the same
   *  values in model and sut */
  static void test_multiply_scalar(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether division by a scalar yields the same
   *  values in model and sut */
  static void test_divide_scalar(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether addition of another arbitrary matrix gives
   *  rise to the same results in model and sut.
   */
  template <typename OtherMatrix>
  static void test_add(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether subtraction of another arbitrary matrix gives
   *  rise to the same results in model and sut.
   */
  template <typename OtherMatrix>
  static void test_subtract(
        const compmat_type& model, const sutmat_type& sut,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default);

  /** Test whether the multiplication with another arbitrary matrix
   *  gives rise to the same result in model and sut.
   */
  template <typename OtherMatrix>
  lazyten_declare_comptest(test_multiply_by);

  /** Test whether the mmult function gives rise to the
   *  expected result for normal application mode
   */
  lazyten_declare_comptest(test_mmult);

  /** Test whether the mmult function gives rise to the
   *  expected result for transposed application mode.
   */
  lazyten_declare_comptest(test_transposed_mmult);

  /** Test whether applying an arbitrary multivector to the matrix
   *  yields the same result in both model and sut */
  template <typename Vector>
  lazyten_declare_comptest(test_apply_to);

  /** Test whether applying an arbitrary multivector to the transpose of the
   *  matrix yields the same result in both model and sut */
  template <typename Vector>
  lazyten_declare_comptest(test_transpose_apply_to);

  /** Test inverse application */
  template <typename Vector>
  lazyten_declare_comptest(test_inv_apply_to);

  /** Are easy problems allowed (by default yes) */
  static bool allow_easy_problems;

  /** Switch to hard problems */
  static void skip_easy_cases() { allow_easy_problems = false; }

 private:
  /** Generator for the extractor ranges in extract_block
   *
   * If numelements==1 we can only return 0
   * Else:
   *   - For easy problems: Try to return more values which are larger than 0
   *   - For hard problems: Always return values larger than 0
   * */
  static Gen<size_type> gen_range_length(size_type numelements) {
    if (numelements == 1) {
      return gen::just<size_type>(0);
    } else {
      if (allow_easy_problems) {
        return gen::weightedOneOf<size_type>(
              {{1, gen::just<size_type>(0)},
               {4, gen::inRange<size_type>(1, numelements)}});
      } else {
        return gen::inRange<size_type>(1, numelements);
      }
    }
  }

  /** Generator for the c_this coefficients for extract_block, mmult and apply
   *
   * Try to return far more more non-zeros than zeros.
   */
  static Gen<scalar_type> gen_c_this() {
    if (allow_easy_problems) {
      return gen::weightedOneOf<scalar_type>(
            {{2, gen::scale(0.9, gen::numeric_around<scalar_type>(1.0))},
             {3, gen::scale(0.9, gen::numeric<scalar_type>())},
             {1, gen::just<scalar_type>(0.0)}});
    } else {
      return gen::scale(0.9, gen::numeric<scalar_type>());
    }
  }

  /** Generator for the c_out coefficients for extract_block, mmult and
   * apply
   *
   * Try to return mostly evenly distributed zeros and non-zeros
   * with focus on some more zeros unless allow_zero is false or we
   * deal with only tough problems, then always return non-zero.
   * */
  static Gen<scalar_type> gen_c_out(const bool allow_zero = true) {
    if (allow_zero && allow_easy_problems) {
      return gen::weightedOneOf<scalar_type>(
            {{2, gen::scale(0.9, gen::numeric_around(1.0))}, {3, gen::just(0.0)}});
    } else {
      return gen::scale(0.9, gen::numeric_around(1.0));
    }
  }

  /** Generator to generate the problem matrices or vectors for
   * apply and mmult
   *
   * If allow_easy_problems is false, we mostly return tensors with more
   * non-zero elements.
   * */
  template <typename Object>
  static Gen<Object> gen_problem_matrix(size_type n_rows, size_type n_cols) {
    if (allow_easy_problems) {
      return gen::numeric_tensor<Object>(n_rows, n_cols);
    } else {
      return gen::scale(4.0, gen::numeric_tensor<Object>(n_rows, n_cols));
    }
  }
};

// Define default value to true
template <typename Model, typename Sut>
bool ComparativeTests<Model, Sut>::allow_easy_problems = true;

//
// -------------------------------------------------------
//

template <typename Matrix1, typename Matrix2, typename Matrix3>
void matrix_product(const Matrix1& lhs, const Matrix2& rhs, Matrix3& out) {
  typedef typename Matrix1::size_type size_type;
  typedef typename Matrix1::scalar_type scalar_type;

  static_assert(
        std::is_same<typename Matrix1::scalar_type, typename Matrix2::scalar_type>::value,
        "All matrices need to have the same scalar type.");
  static_assert(
        std::is_same<typename Matrix2::scalar_type, typename Matrix3::scalar_type>::value,
        "All matrices need to have the same scalar type.");

  static_assert(
        std::is_same<typename Matrix1::size_type, typename Matrix2::size_type>::value,
        "All matrices need to have the same size type.");
  static_assert(
        std::is_same<typename Matrix2::size_type, typename Matrix3::size_type>::value,
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

template <typename Matrix1, typename Matrix2>
void matrix_transpose(const Matrix1& in, Matrix2& trans) {
  typedef typename Matrix1::size_type size_type;

  static_assert(
        std::is_same<typename Matrix1::scalar_type, typename Matrix2::scalar_type>::value,
        "All matrices need to have the same scalar type.");
  static_assert(
        std::is_same<typename Matrix1::size_type, typename Matrix2::size_type>::value,
        "All matrices need to have the same size type.");
  assert_size(in.n_rows(), trans.n_cols());
  assert_size(in.n_cols(), trans.n_rows());

  for (size_type row = 0; row < in.n_rows(); ++row) {
    for (size_type col = 0; col < in.n_cols(); ++col) {
      trans(col, row) = in(row, col);
    }
  }
}

template <typename Matrix, typename Vector1, typename Vector2>
void matrix_apply(const Matrix& lhs, const Vector1& rhs, Vector2& out) {
  typedef typename Matrix::size_type size_type;
  typedef typename Matrix::scalar_type scalar_type;

  static_assert(
        std::is_same<typename Matrix::scalar_type, typename Vector1::scalar_type>::value,
        "All matrices/vectors need to have the same scalar type.");
  static_assert(
        std::is_same<typename Matrix::scalar_type, typename Vector2::scalar_type>::value,
        "All matrices/vectors need to have the same scalar type.");

  static_assert(
        std::is_same<typename Matrix::size_type, typename Vector1::size_type>::value,
        "All matrices/vectors need to have the same size type.");
  static_assert(
        std::is_same<typename Matrix::size_type, typename Vector2::size_type>::value,
        "All matrices/vectors need to have the same size type.");

  assert_size(lhs.n_rows(), out.size());
  assert_size(lhs.n_cols(), rhs.size());

  for (size_type i = 0; i < lhs.n_rows(); ++i) {
    scalar_type res = 0.;
    for (size_type j = 0; j < lhs.n_cols(); ++j) {
      res += lhs(i, j) * rhs(j);
    }
    out(i) = res;
  }
}

//
// ----------
//
//
template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_equivalence(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  // TODO exists in indexable_tests
  //
  // Cannot be used yet, since operator== not overloaded
  RC_ASSERT_NC(model == numcomp(sut).tolerance(tolerance));
}

lazyten_define_comptest(test_copy) {
  // TODO exists in indexable_tests
  idx_cmptests::test_copy(model, sut, tolerance);
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_element_access(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  size_type row = *gen::inRange<size_type>(0, model.n_rows()).as("row index");
  size_type col = *gen::inRange<size_type>(0, model.n_cols()).as("col index");
  RC_ASSERT_NC(model(row, col) == (numcomp(sut(row, col)).tolerance(tolerance)));

  // TODO exists partially in indexable_tests.hh
  idx_cmptests::test_element_access_vectorised(model, sut, tolerance);
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_setting_elements_indexed(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {

  auto modify_row = *gen::inRange<size_type>(0, model.n_rows()).as("Row to modify");
  auto modify_col = *gen::inRange<size_type>(0, model.n_cols()).as("Col to modify");
  auto value = *gen::arbitrary<scalar_type>().as("New value");

  // Keep a copy of the original:
  sutmat_type sut_copy{sut};

  // Modify the value
  sut_copy(modify_row, modify_col) = value;

  // Check it:
  for (auto row : range(model.n_rows())) {
    for (auto col : range(model.n_cols())) {
      if (row == modify_row && col == modify_col) {
        RC_ASSERT_NC(sut_copy(row, col) == numcomp(value).tolerance(tolerance));
      } else {
        RC_ASSERT_NC(sut_copy(row, col) == numcomp(model(row, col)).tolerance(tolerance));
      }
    }
  }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_setting_elements_vectorised(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
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
      RC_ASSERT_NC(sut_copy[i] == numcomp(value).tolerance(tolerance));
    } else {
      RC_ASSERT_NC(sut_copy[i] == numcomp(model[i]).tolerance(tolerance));
    }
  }
}

lazyten_define_comptest(test_extract_block) {
  typedef typename StoredTypeOf<Sut>::type stored_matrix_type;

  size_type nrows =
        *gen_range_length(model.n_rows()).as("Number of rows of the extracted matrix");
  size_type ncols =
        *gen_range_length(model.n_cols()).as("Number of cols of the extracted matrix");
  size_type start_row =
        *gen::inRange<size_type>(0, model.n_rows() - nrows).as("Start row");
  size_type start_col =
        *gen::inRange<size_type>(0, model.n_cols() - ncols).as("Start col");
  auto out = *gen::scale(0.9, gen::numeric_tensor<stored_matrix_type>(nrows, ncols))
                    .as("Matrix to extract values into");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // allow c_M only to be zero if out does not have a small norm.
  auto c_M = *gen_c_out(norm_frobenius(out) >= 1e-12).as("Coefficient for out matrix");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "extract_block: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(c_M == 0.0, "extract_block: Output coefficient is zero(c_M == 0)");
  RC_CLASSIFY(nrows == 0 || ncols == 0, "extract_block: Empty extractor range");
  if (norm_frobenius(out) == 0) {
    RC_CLASSIFY(true, "extract_block: Outmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(out) < 1e-12, "extract_block: Outmatrix has small norm");
  }
#endif

  // Copy the out, that we have an unmodified version
  auto out_copy(out);

  // Extract the block from the matrix
  sut.extract_block(out, start_row, start_col, Transposed::None, c_this, c_M);

  for (size_type i = 0; i < nrows; ++i) {
    for (size_type j = 0; j < ncols; ++j) {
      size_type row = i + start_row;
      size_type col = j + start_col;

      scalar_type model_element = c_this * model(row, col) + c_M * out_copy(i, j);
      RC_ASSERT_NC(out(i, j) == numcomp(model_element).tolerance(tolerance));
    }  // j
  }    // i
}

lazyten_define_comptest(test_extract_transpose_block) {
  typedef typename StoredTypeOf<Sut>::type stored_matrix_type;

  size_type nrows =
        *gen_range_length(model.n_cols()).as("Number of rows of the extracted matrix");
  size_type ncols =
        *gen_range_length(model.n_rows()).as("Number of cols of the extracted matrix");
  size_type start_row = *gen::inRange<size_type>(0, model.n_cols() - nrows)
                               .as("Start row in transp. sut");
  size_type start_col = *gen::inRange<size_type>(0, model.n_rows() - ncols)
                               .as("Start col in transp. sut");
  auto out = *gen::scale(0.9, gen::numeric_tensor<stored_matrix_type>(nrows, ncols))
                    .as("Matrix to extract values into");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // allow c_M only to be zero if out does not have a small norm.
  auto c_M = *gen_c_out(norm_frobenius(out) >= 1e-12).as("Coefficient for out matrix");

// TODO bool conjtransp = *gen::arbitrary<bool>().as("Apply complex
// conjugation");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "extract_tr_block: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(c_M == 0.0, "extract_tr_block: Output coefficient is zero (c_M == 0)");
  RC_CLASSIFY(nrows == 0 || ncols == 0, "extract_tr_block: Empty extractor range");
  if (norm_frobenius(out) == 0) {
    RC_CLASSIFY(true, "extract_tr_block: Outmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(out) < 1e-12,
                "extract_tr_block: Outmatrix has small norm");
  }
#endif

  // Copy the out, that we have an unmodified version
  auto out_copy(out);

  // Extract the block from the matrix
  sut.extract_block(out, start_row, start_col, Transposed::Trans, c_this, c_M);

  for (size_type i = 0; i < nrows; ++i) {
    for (size_type j = 0; j < ncols; ++j) {
      size_type row = i + start_row;
      size_type col = j + start_col;

      scalar_type model_element = c_this * model(col, row) + c_M * out_copy(i, j);
      RC_ASSERT_NC(out(i, j) == numcomp(model_element).tolerance(tolerance));
    }  // j
  }    // i
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_readonly_iterator(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  auto it_const = sut.cbegin();
  auto it = sut.begin();

  for (size_type i = 0; i < model.n_rows(); ++i) {
    for (size_type j = 0; j < model.n_cols(); ++j) {
      // Assert that const and non-const iterators
      // behave identically:
      RC_ASSERT(*it == *it_const);
      RC_ASSERT(it.row() == it_const.row());
      RC_ASSERT(it.col() == it_const.col());

      if (std::abs(model(i, j)) > 10. * std::numeric_limits<scalar_type>::epsilon()) {
        // Larger than zero element -> expect to be present
        RC_ASSERT(it.row() == i);
        RC_ASSERT(it.col() == j);
        RC_ASSERT_NC(model(i, j) == numcomp(*it).tolerance(tolerance));
      }

      if (it.indices() <= std::make_pair(i, j)) {
        // increment iterators if necessary and continue
        ++it;
        ++it_const;
      }
    }
  }

  for (auto it = std::begin(sut); it != std::end(sut); ++it) {
    RC_ASSERT_NC(model(it.row(), it.col()) == numcomp(*it).tolerance(tolerance));
  }
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_readwrite_iterator(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  auto modify_row = *gen::inRange<size_type>(0, model.n_rows()).as("Row to modify");
  auto modify_col = *gen::inRange<size_type>(0, model.n_cols()).as("Col to modify");
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
        RC_ASSERT_NC(sut_copy(i, j) == numcomp(value).tolerance(tolerance));
      } else {
        RC_ASSERT_NC(model(i, j) == numcomp(sut_copy(i, j)).tolerance(tolerance));
      }
    }
  }
}

lazyten_define_comptest(test_norm_l1) {
  scalar_type norm{0};
  for (size_type col = 0; col < model.n_cols(); ++col) {
    scalar_type accu{0};
    for (size_type row = 0; row < model.n_rows(); ++row) {
      accu += std::abs(model(row, col));
    }
    norm = std::max(norm, accu);
  }
  RC_ASSERT_NC(norm_l1(sut) == numcomp(norm).tolerance(tolerance));
}

lazyten_define_comptest(test_norm_linf) {
  scalar_type norm{0};
  for (size_type row = 0; row < model.n_rows(); ++row) {
    scalar_type accu{0};
    for (size_type col = 0; col < model.n_cols(); ++col) {
      accu += std::abs(model(row, col));
    }
    norm = std::max(norm, accu);
  }
  RC_ASSERT_NC(norm_linf(sut) == numcomp(norm).tolerance(tolerance));
}

lazyten_define_comptest(test_norm_frobenius) {
  scalar_type frobenius_squared{0};
  for (size_type i = 0; i < model.n_rows(); ++i) {
    for (size_type j = 0; j < model.n_cols(); ++j) {
      frobenius_squared += model(i, j) * model(i, j);
    }
  }
  RC_ASSERT_NC(norm_frobenius_squared(sut) ==
               numcomp(frobenius_squared).tolerance(tolerance));

  scalar_type frobenius = std::sqrt(frobenius_squared);
  RC_ASSERT_NC(norm_frobenius(sut) == numcomp(frobenius).tolerance(tolerance));
}

lazyten_define_comptest(test_trace) {
  RC_PRE(model.n_rows() == model.n_cols());

  scalar_type res{0};
  for (size_type i = 0; i < model.n_cols(); ++i) {
    res += model(i, i);
  }

  RC_ASSERT_NC(trace(sut) == numcomp(res).tolerance(tolerance));
}

lazyten_define_comptest(test_minmax) {
  // TODO exists in indexable_tests.hh
  idx_cmptests::test_minmax(model, sut, tolerance);
}

lazyten_define_comptest(test_elementwise) {
  // TODO exists in indexable_tests.hh
  idx_cmptests::test_elementwise(model, sut, tolerance);
}

lazyten_define_comptest(test_accumulate) {
  // TODO exists in indexable_tests.hh
  idx_cmptests::test_accumulate(model, sut, tolerance);
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_multiply_scalar(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  // Generate an arbitrary factor, but not too large
  auto c = *gen::numeric<scalar_type>().as("Coefficient");

  // Do the multiplication:
  auto res = sut * c;
  auto res2 = c * sut;

  compmat_type res_model{model.n_rows(), model.n_cols(), false};
  for (auto row : range(model.n_rows())) {
    for (auto col : range(model.n_cols())) {
      res_model(row, col) = model(row, col) * c;
    }
  }

  RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
  RC_ASSERT_NC(res_model == numcomp(res2).tolerance(tolerance));
}

template <typename CompMatrix, typename SutMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_divide_scalar(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  // Generate an arbitrary factor, but not too large
  auto c = *gen::numeric_around<scalar_type>(1.0).as("Coefficient");

  // Do the multiplication:
  auto res = sut / c;

  compmat_type res_model{model.n_rows(), model.n_cols(), false};
  for (auto row : range(model.n_rows())) {
    for (auto col : range(model.n_cols())) {
      res_model(row, col) = model(row, col) / c;
    }
  }
  RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

lazyten_define_comptest_tmpl(test_multiply_by, OtherMatrix) {
  // The size of the other matrix to multiply mlhs with:
  // Note: the RHS of inRange is exclusive
  auto othersize = *gen::numeric_size<2>().as("Number of columns of the rhs matrix");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(othersize == 0, "test_mult_by: Empty rhs matrix");
#endif

  // Generate another matrix:
  OtherMatrix mrhs = *gen::scale(
        0.95,
        gen::fixed_size<OtherMatrix>(model.n_cols(), othersize).as("Multiplication rhs"));

#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(norm_frobenius(mrhs) == 0, "test_mult_by: Zero matrix multiplied");
#endif

  // Do the low-level multiplication:
  Model result_model(model.n_rows(), othersize, false);
  matrix_product(model, mrhs, result_model);

  // multiply them via operator*
  auto result = sut * mrhs;

  // Check that the results are equivalent
  RC_ASSERT_NC(result_model == numcomp(result).tolerance(tolerance));
}

lazyten_define_comptest(test_mmult) {
  typedef typename StoredTypeOf<Sut>::type stored_matrix_type;

  auto col = *gen::numeric_size<2>().as("Number of columns of the rhs matrix");
  auto mrhs =
        *gen::scale(0.95, gen_problem_matrix<stored_matrix_type>(model.n_cols(), col))
               .as("Multiplication matrix");
  auto res =
        *gen::scale(0.85, gen_problem_matrix<stored_matrix_type>(model.n_rows(), col))
               .as("Result matrix");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // allow c_M only to be zero if out does not have a small norm.
  auto coeff = *gen_c_out(norm_frobenius(res) >= 1e-12).as("Coefficient for out matrix");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "mmult: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(coeff == 0.0, "mmult: Output coefficient is zero (c_M == 0)");
  RC_CLASSIFY(col == 0, "mmult: Empty output matrix");

  if (norm_frobenius(mrhs) == 0) {
    RC_CLASSIFY(true, "mmult: Inmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(mrhs) < 1e-12, "mmult: Inmatrix has small norm");
  }
  if (norm_frobenius(res) == 0) {
    RC_CLASSIFY(true, "mmult: Outmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(res) < 1e-12, "mmult: Outmatrix has small norm");
  }
#endif

  // Perform operations on model:
  Model model_mmult(model.n_rows(), col, false);
  Model model_res(res.n_rows(), res.n_cols(), false);
  matrix_product(model, mrhs, model_mmult);
  for (size_type i = 0; i < res.n_rows(); ++i) {
    for (size_type j = 0; j < res.n_cols(); ++j) {
      model_res(i, j) = coeff * res(i, j) + c_this * model_mmult(i, j);
    }
  }
  sut.mmult(mrhs, res, Transposed::None, c_this, coeff);

  RC_ASSERT_NC(res == numcomp(model_res).tolerance(tolerance));
}

lazyten_define_comptest(test_transposed_mmult) {
  typedef typename StoredTypeOf<Sut>::type stored_matrix_type;

  auto col = *gen::numeric_size<2>().as("Number of columns of the rhs matrix");
  auto mrhs =
        *gen::scale(0.95, gen_problem_matrix<stored_matrix_type>(model.n_rows(), col))
               .as("Multiplication matrix");
  auto res =
        *gen::scale(0.85, gen_problem_matrix<stored_matrix_type>(model.n_cols(), col))
               .as("Result matrix");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // allow c_M only to be zero if out does not have a small norm.
  auto coeff = *gen_c_out(norm_frobenius(res) >= 1e-12).as("Coefficient for out matrix");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "mmult_trans: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(coeff == 0.0, "mmult_trans: Output coefficient is zero (c_M == 0)");
  RC_CLASSIFY(col == 0, "mmult_trans: Empty output matrix");
  if (norm_frobenius(mrhs) == 0) {
    RC_CLASSIFY(true, "mmult_trans: Inmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(mrhs) < 1e-12, "mmult_trans: Inmatrix has small norm");
  }
  if (norm_frobenius(res) == 0) {
    RC_CLASSIFY(true, "mmult_trans: Outmatrix is zero");
  } else {
    RC_CLASSIFY(norm_frobenius(res) < 1e-12, "mmult_trans: Outmatrix is zero");
  }
#endif

  // Perform operations on model:
  Model model_transposed(model.n_cols(), model.n_rows());
  Model model_mmult(model.n_cols(), col, false);
  Model model_res(res.n_rows(), res.n_cols(), false);
  matrix_transpose(model, model_transposed);
  matrix_product(model_transposed, mrhs, model_mmult);
  for (size_type i = 0; i < res.n_rows(); ++i) {
    for (size_type j = 0; j < res.n_cols(); ++j) {
      model_res(i, j) = coeff * res(i, j) + c_this * model_mmult(i, j);
    }
  }
  sut.mmult(mrhs, res, Transposed::Trans, c_this, coeff);

  RC_ASSERT_NC(res == numcomp(model_res).tolerance(tolerance));
}

lazyten_define_comptest_tmpl(test_apply_to, Vector) {
  auto n_vectors = *gen::inRange<size_type>(1, 6).as("Number of vectors to apply to");
  auto mvin = *gen_problem_matrix<MultiVector<Vector>>(model.n_cols(), n_vectors)
                     .as("Vectors to apply to");
  auto mvout = *gen::scale(0.9, gen_problem_matrix<MultiVector<Vector>>(model.n_rows(),
                                                                        n_vectors))
                      .as("Result vectors");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // Norm of outvectors:
  const auto l2norms_out = norm_l2_squared(mvout);
  const scalar_type norm_out =
        std::sqrt(std::accumulate(std::begin(l2norms_out), std::end(l2norms_out), 0));

  // Norm of invectors:
  const auto l2norms_in = norm_l2_squared(mvin);
  const scalar_type norm_in =
        std::sqrt(std::accumulate(std::begin(l2norms_in), std::end(l2norms_in), 0));

  auto coeff = *gen_c_out(norm_out >= 1e-12).as("Coefficient for out vectors");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "apply: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(coeff == 0.0, "apply: Output coefficient is zero (c_y == 0)");
  if (norm_in == 0) {
    RC_CLASSIFY(true, "apply: Invectors are zero ");
  } else {
    RC_CLASSIFY(norm_in < 1e-12, "apply: Invectors have small norm");
  }
  if (norm_out == 0) {
    RC_CLASSIFY(true, "apply: Outvectors are zero");
  } else {
    RC_CLASSIFY(norm_out < 1e-12, "apply: Outvectors have small norm");
  }
#endif

  // Perform operation on model:
  auto model_out = mvout.copy_deep();
  for (size_type vec = 0; vec < n_vectors; ++vec) {
    const auto& vi = mvin[vec];
    const auto& vo = mvout[vec];
    auto& mo = model_out[vec];

    Vector tmp(vo.size());
    matrix_apply(model, vi, tmp);
    mo = c_this * tmp + coeff * vo;
  }

  sut.apply(mvin, mvout, Transposed::None, c_this, coeff);

  RC_ASSERT_NC(model_out == numcomp(mvout).tolerance(tolerance));
}

lazyten_define_comptest_tmpl(test_transpose_apply_to, Vector) {
  auto n_vectors = *gen::inRange<size_type>(1, 6).as("Number of vectors to apply to");
  auto mvin = *gen_problem_matrix<MultiVector<Vector>>(model.n_rows(), n_vectors)
                     .as("Vectors to apply to");
  auto mvout = *gen::scale(0.9, gen_problem_matrix<MultiVector<Vector>>(model.n_cols(),
                                                                        n_vectors))
                      .as("Result vectors");
  auto c_this = *gen_c_this().as("Coefficient for sut matrix");

  // Norm of outvectors:
  const auto l2norms_out = norm_l2_squared(mvout);
  const scalar_type norm_out =
        std::sqrt(std::accumulate(std::begin(l2norms_out), std::end(l2norms_out), 0));

  // Norm of invectors:
  const auto l2norms_in = norm_l2_squared(mvin);
  const scalar_type norm_in =
        std::sqrt(std::accumulate(std::begin(l2norms_in), std::end(l2norms_in), 0));

  auto coeff = *gen_c_out(norm_out >= 1e-12).as("Coefficient for out vectors");
#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(c_this == 0.0, "apply_trans: Matrix coefficient is zero (c_this==0)");
  RC_CLASSIFY(coeff == 0.0, "apply_trans: Output coefficient is zero (c_y == 0)");
  if (norm_in == 0) {
    RC_CLASSIFY(true, "apply_trans: Invectors are zero");
  } else {
    RC_CLASSIFY(norm_in < 1e-12, "apply_trans: Invectors have small norm");
  }
  if (norm_out == 0) {
    RC_CLASSIFY(true, "apply_trans: Outvectors are zero");
  } else {
    RC_CLASSIFY(norm_out < 1e-12, "apply_trans: Outvectors have small norm");
  }
#endif

  // Perform operation on model:
  Model modeltransp(model.n_cols(), model.n_rows());
  matrix_transpose(model, modeltransp);
  auto model_out = mvout.copy_deep();
  for (size_type vec = 0; vec < n_vectors; ++vec) {
    const auto& vi = mvin[vec];
    const auto& vo = mvout[vec];
    auto& mo = model_out[vec];

    Vector tmp(vo.size());
    matrix_apply(modeltransp, vi, tmp);
    mo = c_this * tmp + coeff * vo;
  }

  sut.apply(mvin, mvout, Transposed::Trans, c_this, coeff);
  RC_ASSERT_NC(model_out == numcomp(mvout).tolerance(tolerance));
}

lazyten_define_comptest_tmpl(test_inv_apply_to, Vector) {
  auto n_vectors = *gen::inRange<size_type>(1, 6).as("Number of vectors to apply to");
  auto ivec_gen =
        gen::with_l2_norm_in_range(0, 1e2, gen::numeric_tensor<Vector>(model.n_cols()));
  auto mvin = *gen::numeric_tensor<MultiVector<Vector>>(n_vectors, ivec_gen)
                     .as("Input vectors");
  auto ovec_gen =
        gen::with_l2_norm_in_range(0, 1e2, gen::numeric_tensor<Vector>(model.n_cols()));
  auto mvout = *gen::numeric_tensor<MultiVector<Vector>>(n_vectors, ovec_gen)
                      .as("Result vectors");
  auto c_this = *gen::map(gen::inRange<long>(10000000, 1000000001), [](long l) {
                   return l / 1e8;
                 }).as("Coefficient for sut matrix");
  auto coeff = *gen_c_out().as("Coefficient for out vectors");

  // Vector for the results:
  MultiVector<Vector> res(mvin.n_elem(), mvin.n_vectors(), false);

  // TODO also test ConjTrans and Trans (should work with this test setup as well!)
  const Transposed transop = Transposed::None;

  auto mvoutcopy = mvout.copy_deep();

  // First apply_inverse:
  try {
    sut.apply_inverse(mvin, mvout, transop, c_this, coeff);
  } catch (const krims::ExcDisabled& e) {
    // This means that this matrix does not support values of coeff other than 0
    coeff = 0;
    sut.apply_inverse(mvin, mvout, transop, c_this, coeff);
  }
  // mvout now holds:
  //   c_this * M^{-1} * mvin + coeff * mvout

  // Now undo that:
  sut.apply(mvout, res, transop, 1 / c_this, 0);

  // Res now holds:
  //  1/c_this * M * mvout = 1/c_this * M * (c_this * M^{-1} * mvin + coeff * mvout)
  //                       = mvin + coeff / c_this * M * mvout

  // Use the model to compute the part we need to add to mvin
  MultiVector<Vector> model_out(mvin.n_elem(), mvin.n_vectors(), false);
  for (size_type vec = 0; vec < n_vectors; ++vec) {
    const auto& vi = mvin[vec];
    const auto& vo = mvoutcopy[vec];
    auto& mo = model_out[vec];

    Vector tmp(mo.size());
    matrix_apply(model, vo, tmp);
    mo = vi + coeff / c_this * tmp;
  }

  // Check that it agrees:
  RC_ASSERT_NC(model_out == numcomp(res).tolerance(tolerance));
}

/** Test whether addition of another arbitrary matrix
 * gives rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_add(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  // generate another matrix of the same size:
  auto madd =
        *gen::fixed_size<OtherMatrix>(model.n_rows(), model.n_cols()).as("Matrix to add");

#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(norm_frobenius(madd) == 0, "add: Zero matrix added");
#endif

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
  RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

/** Test whether subtraction of another arbitrary matrix
 * gives
 *  rise to the same results in model and sut.
 */
template <typename CompMatrix, typename SutMatrix>
template <typename OtherMatrix>
void ComparativeTests<CompMatrix, SutMatrix>::test_subtract(
      const compmat_type& model, const sutmat_type& sut,
      const NumCompAccuracyLevel tolerance) {
  // generate another matrix of the same size:
  auto msub = *gen::fixed_size<OtherMatrix>(model.n_rows(), model.n_cols())
                     .as("Matrix to subtract");

#ifdef LAZYTEN_TESTS_VERBOSE
  RC_CLASSIFY(norm_frobenius(msub) == 0, "sub: Zero matrix subtracted");
#endif

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
  RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

}  // matrix_tests
}  // namespace tests
}  // namespace lazyten
