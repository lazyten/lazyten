//
// Copyright (C) 2016-17 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

#include "lazy_matrix_tests.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <linalgwrap/LazyMatrix_i.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/TestingUtils.hh>
#include <linalgwrap/inverse.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;
using namespace krims;

namespace inverse_tests {
template <typename Matrix>
class SimpleInvertible : public LazyMatrix_i<Matrix> {
 public:
  typedef typename Matrix::scalar_type scalar_type;
  typedef typename LazyMatrix_i<Matrix>::lazy_matrix_expression_ptr_type
        lazy_matrix_expression_ptr_type;

  void update(const krims::GenMap& /*map*/) override {}

  bool has_apply_inverse() const override { return true; }

  size_t n_rows() const override { return m_inner.n_rows(); }
  size_t n_cols() const override { return m_inner.n_cols(); }
  scalar_type operator()(size_t row, size_t col) const override {
    return m_inner(row, col);
  }

  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<SimpleInvertible, VectorIn, VectorOut>...>
  void apply_inverse(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1, const scalar_type c_y = 0) const {
    m_inv.apply(x, y, mode, c_this, c_y);
  }

  void apply_inverse(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
                     MultiVector<MutableMemoryVector_i<scalar_type>>& y,
                     const Transposed mode = Transposed::None,
                     const scalar_type c_this = 1,
                     const scalar_type c_y = 0) const override {
    m_inv.apply(x, y, mode, c_this, c_y);
  }

  template <typename VectorIn, typename VectorOut,
            mat_vec_apply_enabled_t<SimpleInvertible, VectorIn, VectorOut>...>
  void apply(const MultiVector<VectorIn>& x, MultiVector<VectorOut>& y,
             const Transposed mode,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const {
    m_inner.apply(x, y, mode, c_this, c_y);
  }

  void apply(const MultiVector<const MutableMemoryVector_i<scalar_type>>& x,
             MultiVector<MutableMemoryVector_i<scalar_type>>& y,
             const Transposed mode = Transposed::None,
             const scalar_type c_this = Constants<scalar_type>::one,
             const scalar_type c_y = Constants<scalar_type>::zero) const override {
    m_inner.apply(x, y, mode, c_this, c_y);
  }

  /** \brief Clone the expression */
  lazy_matrix_expression_ptr_type clone() const override {
    // return a copy enwrapped in the pointer type
    return lazy_matrix_expression_ptr_type(new SimpleInvertible(*this));
  }

 private:
  Matrix m_inner{
        {-20.916781826309034, -44.00752036167732, -6.51669454312858, 11.238980631660937,
         -18.735649841069375, -4.28414790510692, -5.979058937790342, -67.15332851719711,
         12.40302313915069, -29.70975657144828},
        {-44.00752036167732, 43.33168132725473, 36.15697968315908, 3.7913373468127958,
         1.067603854903723, -32.683319038449795, 20.538368062868216, 89.8009658926145,
         -36.059439706338054, -11.420464419731644},
        {-6.51669454312858, 36.15697968315908, 36.823522397823695, 12.385902439431447,
         62.87036215328271, -11.480111108312158, -29.070550765582993, -18.660443732773743,
         -26.588676677121526, -15.82246819798398},
        {11.238980631660937, 3.7913373468127958, 12.385902439431447, -47.71687996968947,
         66.94552735606189, -53.624587170585926, 44.68377336085291, 27.954146861225183,
         50.84098704687406, -33.6574104715489},
        {-18.735649841069375, 1.067603854903723, 62.87036215328271, 66.94552735606189,
         42.51698652146962, 20.073997091362173, -4.887610576472497, 61.09603368533391,
         38.96748290019141, -5.6152672035205455},
        {-4.28414790510692, -32.683319038449795, -11.480111108312158, -53.624587170585926,
         20.073997091362173, 56.07657950843145, -6.888794224294372, -36.093302839619525,
         49.957587316469386, -78.61170150098152},
        {-5.979058937790342, 20.538368062868216, -29.070550765582993, 44.68377336085291,
         -4.887610576472497, -6.888794224294372, -60.84292609644129, -40.90038300715179,
         8.024898569242339, 22.116879619681868},
        {-67.15332851719711, 89.8009658926145, -18.660443732773743, 27.954146861225183,
         61.09603368533391, -36.093302839619525, -40.90038300715179, 22.25841518743931,
         56.36728387202521, -7.83020507279663},
        {12.40302313915069, -36.059439706338054, -26.588676677121526, 50.84098704687406,
         38.96748290019141, 49.957587316469386, 8.024898569242339, 56.36728387202521,
         20.320675131824714, -4.494980491142471},
        {-29.70975657144828, -11.420464419731644, -15.82246819798398, -33.6574104715489,
         -5.6152672035205455, -78.61170150098152, 22.116879619681868, -7.83020507279663,
         -4.494980491142471, 51.908640287424475}};

  Matrix m_inv{{-0.00501497646750753, -0.0003100248081494239, -0.0024005720015216965,
                0.01140669034137758, -0.0015238504900875935, -0.0037934390620201464,
                0.0123174608770872, -0.007875486044424788, -0.0018822351788730575,
                -0.008783013339318657},
               {-0.0003100248081494201, -0.007064218837994684, -0.0013017191398094002,
                -0.0020871850097904103, -0.004659232540056946, -0.009202725632069126,
                -0.01409338454832073, 0.011348296299031526, -0.005545441385712342,
                -0.010686145906823888},
               {-0.0024005720015217, -0.0013017191398094, 0.0020234275306485025,
                -0.0009757717848136165, 0.010649791959989274, 0.0018357665575258742,
                -0.0009797072653025853, -0.003987713226331778, -0.00895621163616647,
                0.0012962512189738217},
               {0.011406690341377585, -0.0020871850097904108, -0.0009757717848136158,
                0.0015608016653788008, -0.0008556615578333949, -0.011756840154112533,
                -0.004069591185721513, 0.003926937337172913, 0.0035350338937473343,
                -0.008481013714945067},
               {-0.00152385049008759, -0.004659232540056942, 0.010649791959989274,
                -0.0008556615578333933, -0.0012977886707000812, 0.0006179862993093697,
                -0.005876602442439707, 0.002247261200966867, 0.006907036888666075,
                0.005030606842483763},
               {-0.0037934390620201447, -0.009202725632069121, 0.001835766557525871,
                -0.011756840154112531, 0.000617986299309369, 0.0010894081166700508,
                -0.017349663484493186, 0.007015192926469238, 0.0011103821842738662,
                -0.0009961278179763675},
               {0.01231746087708721, -0.014093384548320724, -0.0009797072653025887,
                -0.004069591185721512, -0.005876602442439715, -0.017349663484493197,
                -0.03529505120488893, 0.014738344834724871, 0.0014689514967670316,
                -0.008509936084101444},
               {-0.007875486044424791, 0.01134829629903152, -0.003987713226331778,
                0.003926937337172912, 0.0022472612009668685, 0.007015192926469239,
                0.01473834483472487, -0.007138124944921143, 0.003009824948191947,
                0.0030912823818593105},
               {-0.0018822351788730612, -0.005545441385712341, -0.00895621163616647,
                0.0035350338937473326, 0.006907036888666078, 0.0011103821842738675,
                0.0014689514967670307, 0.003009824948191948, -0.005152375903511174,
                -0.0009244813044098285},
               {-0.008783013339318661, -0.010686145906823884, 0.0012962512189738187,
                -0.008481013714945067, 0.005030606842483765, -0.000996127817976364,
                -0.008509936084101437, 0.0030912823818593097, -0.0009244813044098311,
                0.009830405383244515}};
};

}  // namespace inverse_tests

TEST_CASE("inverse function", "[inverse]") {
  using namespace inverse_tests;

  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> matrix_type;
  typedef matrix_type::vector_type vector_type;
  const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Sloppy;

  SECTION("Test with SimpleInvertible") {
    auto testable = [] {
      SimpleInvertible<matrix_type> m;
      auto m_inv = inverse(m);

      auto op = m * m_inv;
      auto op2 = m_inv * m;

      auto vec_gen = gen::with_l2_norm_in_range(
            0, 1e4, gen::numeric_tensor<vector_type>(m.n_cols()));
      auto n_vecs =
            *gen::distinctFrom(gen::numeric_size<2>().as("Number of vectors"), 0ul);
      auto mv =
            *gen::numeric_tensor<MultiVector<vector_type>>(n_vecs, vec_gen).as("vectors");

      auto res = op * mv;
      auto res2 = op2 * mv;

      RC_ASSERT_NC(res == numcomp(mv).tolerance(tolerance));
      RC_ASSERT_NC(res2 == numcomp(mv).tolerance(tolerance));
    };
    REQUIRE(rc::check("Test with SimpleInvertable", testable));
  }  // inverse * matrix == identity

  auto random_matrix_generator = [] {
    auto n = *gen::scale(0.8, gen::numeric_size<2>()).as("Matrix size");
    auto mat_gen = gen::with_properties(gen::numeric_tensor<matrix_type>(n, n),
                                        OperatorProperties::RealSymmetric);
    auto add_one = [](matrix_type m) {
      // Add one on the diagonal to make matrix less singular:
      for (size_t i = 0; i < m.n_rows(); ++i) {
        m(i, i) += 1.;
      }
      return m;
    };
    auto mat = *gen::map(mat_gen, add_one).as("Problem matrix");

    return mat;
  };

  SECTION("Test with arbitrary hermitian matrix and make_invertible") {
    auto testable = [random_matrix_generator]() {
      auto m = random_matrix_generator();
      auto m_inv = make_invertible(m);

      auto op = m_inv * inverse(m_inv);
      auto op2 = inverse(m_inv) * m_inv;

      auto vec_gen = gen::with_l2_norm_in_range(
            0, 1e4, gen::numeric_tensor<vector_type>(m_inv.n_cols()));
      auto n_vecs =
            *gen::distinctFrom(gen::numeric_size<2>().as("Number of vectors"), 0ul);
      auto mv =
            *gen::numeric_tensor<MultiVector<vector_type>>(n_vecs, vec_gen).as("vectors");

      auto res = op * mv;
      auto res2 = op2 * mv;

      RC_ASSERT_NC(res == numcomp(mv).tolerance(tolerance));
      RC_ASSERT_NC(res2 == numcomp(mv).tolerance(tolerance));
    };

    REQUIRE(rc::check("Test with arbitrary hermitian matrix and make_invertible",
                      testable));
  }

  SECTION("Test object resulting from make_invertible.") {
    // Generator for the model.
    auto model_generator = [](matrix_type m) { return m; };

    // Generator for the sut
    struct InvertibleGenerator {
      auto operator()(matrix_type m) -> decltype(make_invertible(m)) {
        m_ptr.reset(new matrix_type(std::move(m)));
        return make_invertible(*m_ptr);
      }
      std::shared_ptr<matrix_type> m_ptr;
    };

    typedef lazy_matrix_tests::TestingLibrary<
          typename std::result_of<InvertibleGenerator(matrix_type)>::type, matrix_type>
          testlib;

    // Make testing a little easier for this section
    auto highertol = NumCompConstants::change_temporary(
          10. * krims::NumCompConstants::default_tolerance_factor);

    testlib tl{random_matrix_generator, model_generator, InvertibleGenerator{},
               "make_invertible(): "};
    tl.enable_inverse_apply();
    tl.run_checks();
  }

}  // namespace inverse
}  // namespace tests
}  // namespace linalgwrap
