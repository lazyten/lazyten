//
// Copyright (C) 2016 by the linalgwrap authors
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

#pragma once
#include "vector_tests.hh"

namespace linalgwrap {
namespace tests {
namespace mutable_vector_tests {
using namespace rc;
using namespace krims;

/** \brief Standard test functions which test a certain
 *  functionality by executing it in the Sut indexable
 *  and in a Model indexable and comparing the results
 *  afterwards.
 *
 *  Tests for mutable vectors.
 *
 *  \tparam Model  Model indexable used for comparison
 *  \tparam Sut   System under test indexable.
 *                      The thing we test.
 **/
template <typename Model, typename Sut>
struct ComparativeTests : public vector_tests::ComparativeTests<Model, Sut> {
    typedef vector_tests::ComparativeTests<Model, Sut> base_type;
    typedef typename base_type::model_type model_type;
    typedef typename base_type::sut_type sut_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::real_type real_type;

    /** Test read-write element access via () at a random place.
     *  compare result against a model at all places except the
     *  changed entry.
     */
    linalgwrap_declare_comptest(test_setting_elements);

    /** Test read-write element access via [] at a random place.
     *  compare result against a model at all places except the
     *  changed entry.
     */
    linalgwrap_declare_comptest(test_setting_elements_vectorised);

    /** Test whether modifying element entries using an iterator
     *  works.
     **/
    linalgwrap_declare_comptest(test_readwrite_iterator);

    /** Test whether multiplication by a scalar yields the same
     *  values in model and sut */
    linalgwrap_declare_comptest(test_multiply_scalar);

    /** Test whether division by a scalar yields the same
     *  values in model and sut */
    linalgwrap_declare_comptest(test_divide_scalar);

    /** Test whether addition of another arbitrary vector gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherVector>
    linalgwrap_declare_comptest(test_add);

    /** Test whether subtraction of another arbitrary vector gives
     *  rise to the same results in model and sut.
     */
    template <typename OtherVector>
    linalgwrap_declare_comptest(test_subtract);

    /** Run all comparative tests (including the ones from vector_tests and
     * indexable_tests) */
    template <typename Args>
    static void run_all(const RCTestableGenerator<Model, Sut, Args>& gen,
                        const std::string& prefix);
};

//
// -------------------------------------------------------
//

linalgwrap_define_comptest(test_setting_elements) {
    auto modify_index =
          *gen::inRange<size_type>(0, model.n_elem()).as("Element to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    sut_type sut_copy{sut};

    // Modify the value
    sut_copy(modify_index) = value;

    // Check it:
    for (size_type i = 0; i < model.n_elem(); ++i) {
        if (i == modify_index) {
            RC_ASSERT_NC(sut_copy(i) == numcomp(value).tolerance(tolerance));
        } else {
            RC_ASSERT_NC(sut_copy(i) == numcomp(model(i)).tolerance(tolerance));
        }
    }
}

linalgwrap_define_comptest(test_setting_elements_vectorised) {
    auto modify_index =
          *gen::inRange<size_type>(0, model.n_elem()).as("Element to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    sut_type sut_copy{sut};

    // Modify the value
    sut_copy[modify_index] = value;

    // Check it:
    for (size_type i = 0; i < model.n_elem(); ++i) {
        if (i == modify_index) {
            RC_ASSERT_NC(sut_copy[i] == numcomp(value).tolerance(tolerance));
        } else {
            RC_ASSERT_NC(sut_copy[i] == numcomp(model[i]).tolerance(tolerance));
        }
    }
}

linalgwrap_define_comptest(test_readwrite_iterator) {
    auto modify_index =
          *gen::inRange<size_type>(0, model.n_elem()).as("Element to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    sut_type sut_copy{sut};

    auto it = sut_copy.begin();
    std::advance(it, modify_index);
    *it = value;

    for (size_type i = 0; i < model.n_elem(); ++i) {
        if (i == modify_index) {
            RC_ASSERT_NC(sut_copy[i] == numcomp(value).tolerance(tolerance));
        } else {
            RC_ASSERT_NC(sut_copy[i] == numcomp(model[i]).tolerance(tolerance));
        }
    }
}

linalgwrap_define_comptest(test_multiply_scalar) {
    // Generate an arbitrary factor, but not too large
    auto c = *gen::numeric<scalar_type>().as("Coefficient");

    // Do the multiplication:
    auto res = sut * c;
    auto res2 = c * sut;

    model_type res_model(model.n_elem(), false);
    for (size_type i = 0; i < model.n_elem(); ++i) {
        res_model[i] = c * model[i];
    }

    RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
    RC_ASSERT_NC(res_model == numcomp(res2).tolerance(tolerance));
}

linalgwrap_define_comptest(test_divide_scalar) {
    // Generate an arbitrary factor, but not too large
    auto c = *gen::numeric_nonZero<scalar_type>().as("Coefficient");

    // Do the multiplication:
    auto res = sut / c;

    model_type res_model(model.n_elem(), false);
    for (size_type i = 0; i < model.n_elem(); ++i) {
        res_model[i] = model[i] / c;
    }

    RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

linalgwrap_define_comptest_tmpl(test_add, OtherVector) {
    auto madd =
          *gen::numeric_tensor<OtherVector>(model.size()).as("Vector to add");

    // Perform the operation
    auto res = sut + madd;

    // and on the model:
    model_type res_model(model.n_elem(), false);
    for (size_type i = 0; i < model.n_elem(); ++i) {
        res_model[i] = model[i] + madd[i];
    }

    // Check that the results are equivalent
    RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

linalgwrap_define_comptest_tmpl(test_subtract, OtherVector) {
    auto msub = *gen::numeric_tensor<OtherVector>(model.size())
                       .as("Vector to substract");

    // Perform the operation
    auto res = sut - msub;

    // and on the model:
    model_type res_model(model.n_elem(), false);
    for (size_type i = 0; i < model.n_elem(); ++i) {
        res_model[i] = model[i] - msub[i];
    }

    // Check that the results are equivalent
    RC_ASSERT_NC(res_model == numcomp(res).tolerance(tolerance));
}

template <typename Model, typename Sut>
template <typename Args>
void ComparativeTests<Model, Sut>::run_all(
      const RCTestableGenerator<Model, Sut, Args>& gen,
      const std::string& prefix) {
    base_type::run_all(gen, prefix);

    const NumCompAccuracyLevel deflvl = NumCompAccuracyLevel::Default;
    const NumCompAccuracyLevel sloppy = NumCompAccuracyLevel::Sloppy;
    const NumCompAccuracyLevel eps = NumCompAccuracyLevel::MachinePrecision;

    // Read-write element access
    CHECK(gen.run_test(prefix + "Altering elements via ()",
                       test_setting_elements, eps));
    CHECK(gen.run_test(prefix + "Altering elements via []",
                       test_setting_elements_vectorised, eps));
    CHECK(gen.run_test(prefix + "Altering elements via iterator",
                       test_readwrite_iterator, eps));

    // Operations
    CHECK(gen.run_test(prefix + "Multiplication by scalar",
                       test_multiply_scalar));
    CHECK(gen.run_test(prefix + "Divide by scalar", test_divide_scalar,
                       base_type::cplx ? sloppy : deflvl));
    CHECK(gen.run_test(prefix + "Add a vector", test_add<Sut>));
    CHECK(gen.run_test(prefix + "Subtract a vector", test_subtract<Sut>));
}

}  // namespace mutable_vector_tests
}  // namespace tests
}  // namespace linalgwrap
