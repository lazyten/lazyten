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
#include "indexable_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;
using namespace krims;

/** Namespace for the components for standard vector tests.
 **/
namespace vector_tests {

/** \brief Standard test functions which test a certain
 *  functionality by executing it in the SutVector
 *  and in a ModelVector and comparing the results
 *  afterwards.
 *
 *  \tparam Model  Model vector used for comparison
 *  \tparam Sut   System under test vector.
 *                      The thing we test.
 **/
template <typename Model, typename Sut>
struct ComparativeTests : public indexable_tests::ComparativeTests<Model, Sut> {
    typedef indexable_tests::ComparativeTests<Model, Sut> base_type;
    typedef typename base_type::model_type model_type;
    typedef typename base_type::sut_type sut_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::real_type real_type;

    // TODO test swap function!

    /** Test read-only element access via () at random places */
    linalgwrap_declare_comptest(test_element_access);

    /** Test the l1 norm function */
    linalgwrap_declare_comptest(test_norm_l1);

    /** Test the linf norm function */
    linalgwrap_declare_comptest(test_norm_linf);

    /** Test the frobenius norm functions (norm_frobenius and
     * norm_frobenius_squared) */
    linalgwrap_declare_comptest(test_norm_l2);

    /** Run all comparative tests (including the ones from indexable_tests) */
    template <typename Args>
    static void run_all(const RCTestableGenerator<Model, Sut, Args>& gen,
                        const std::string& prefix);
};

//
// -------------------------------------------------------
//

linalgwrap_define_comptest(test_element_access) {
    size_type i = *gen::inRange<size_type>(0u, model.size()).as("index");
    RC_ASSERT_NC(model(i) == numcomp(sut(i)).tolerance(tolerance));
}

linalgwrap_define_comptest(test_norm_l1) {
    real_type norm{0};
    for (size_type i = 0; i < model.n_elem(); ++i) {
        norm += std::abs(model[i]);
    }
    RC_ASSERT_NC(norm_l1(sut) == numcomp(norm).tolerance(tolerance));
}

linalgwrap_define_comptest(test_norm_linf) {
    real_type norm{0};
    for (size_type i = 0; i < model.n_elem(); ++i) {
        norm = std::max(norm, std::abs(model[i]));
    }
    RC_ASSERT_NC(norm_linf(sut) == numcomp(norm).tolerance(tolerance));
}

linalgwrap_define_comptest(test_norm_l2) {
    real_type norm_squared{0};
    for (size_type i = 0; i < model.n_elem(); ++i) {
        norm_squared += std::abs(model[i]) * std::abs(model[i]);
    }
    real_type norm = std::sqrt(norm_squared);

    RC_ASSERT_NC(norm_l2_squared(sut) ==
                 numcomp(norm_squared).tolerance(tolerance));
    RC_ASSERT_NC(norm_l2(sut) == numcomp(norm).tolerance(tolerance));
}

template <typename Model, typename Sut>
template <typename Args>
void ComparativeTests<Model, Sut>::run_all(
      const RCTestableGenerator<Model, Sut, Args>& gen,
      const std::string& prefix) {
    base_type::run_all(gen, prefix);

    const NumCompAccuracyLevel eps = NumCompAccuracyLevel::MachinePrecision;
    const NumCompAccuracyLevel eps10 =
          NumCompAccuracyLevel::TenMachinePrecision;

    CHECK(gen.run_test(prefix + "Element access via ()", test_element_access,
                       eps));

    CHECK(gen.run_test(prefix + "l1 norm", test_norm_l1, eps10));
    CHECK(gen.run_test(prefix + "l2 norm", test_norm_l2, eps10));
    CHECK(gen.run_test(prefix + "linf norm", test_norm_linf, eps10));
}

}  // vector_tests
}  // namespace tests
}  // namespace linalgwrap
