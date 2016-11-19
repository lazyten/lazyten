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
#include "mutable_vector_tests.hh"

namespace linalgwrap {
namespace tests {

/** Namespace for default tests for stored vectors */
namespace stored_vector_tests {

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
struct ComparativeTests
      : public mutable_vector_tests::ComparativeTests<Model, Sut> {
    typedef indexable_tests::ComparativeTests<Model, Sut> base_type;
    typedef typename base_type::model_type model_type;
    typedef typename base_type::sut_type sut_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::real_type real_type;

    /** Run all comparative tests (including the ones from indexable_tests,
     * vector_tests and mutable_vector_tests) */
    template <typename Args>
    static void run_all(const RCTestableGenerator<Model, Sut, Args>& gen,
                        const std::string& prefix) {

        once_test_conversion_from_indexable();
        once_test_initialiser_list_constructor();
        once_test_as_scalar();

        base_type::run_all(gen, prefix);
    }

    static void once_test_initialiser_list_constructor();
    static void once_test_conversion_from_indexable();
    static void once_test_as_scalar();
};

template <typename Model, typename Sut, typename Args>
void run_with_generator(const RCTestableGenerator<Model, Sut, Args>& gen,
                        const std::string prefix = "") {
    ComparativeTests<Model, Sut>::run_all(gen, prefix);
}

//
// ---------------------------------------
//

template <typename Model, typename Sut>
void ComparativeTests<Model, Sut>::once_test_initialiser_list_constructor() {
    sut_type v{11.0, 22.0, 33.0, 44.0};

    CHECK((v.n_elem() == 4));
    CHECK((v.size() == 4));

    for (size_type i = 0; i < v.size(); ++i) {
        CHECK((v(i) == 11. * (i + 1)));
        CHECK((v[i] == 11. * (i + 1)));
    }

    Sut vref(4);
    vref[0] = 11.0;
    vref[1] = 22.0;
    vref[2] = 33.0;
    vref[3] = 44.0;
    CHECK(v == vref);
}

template <typename Model, typename Sut>
void ComparativeTests<Model, Sut>::once_test_conversion_from_indexable() {
    model_type m({1., 2., 3., 4., 5.});
    Sut v(m);

    for (size_type i = 0; i < v.size(); ++i) {
        CHECK(v[i] == m[i]);
    }
}

template <typename Model, typename Sut>
void ComparativeTests<Model, Sut>::once_test_as_scalar() {
    Sut v{42.};
    scalar_type s = as_scalar(v);
    CHECK(s == v[0]);
}

}  // namespace stored_matrix_tests
}  // namespace tests
}  // namespace linalgwrap
