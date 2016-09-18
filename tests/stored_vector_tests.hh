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
#include "RCTestableGenerator.hh"
#include "generators.hh"
#include "vector_tests.hh"
#include <catch.hpp>

namespace linalgwrap {
namespace tests {

/** Namespace for default tests for stored vectors */
namespace stored_vector_tests {

template <typename Scalar>
struct VectorModel : private std::vector<Scalar>, public Vector_i<Scalar> {
    typedef Vector_i<Scalar> base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename std::vector<Scalar> container_type;

    // Bring in important stuff from std::vector
    using typename container_type::iterator;
    using typename container_type::const_iterator;
    using container_type::begin;
    using container_type::cbegin;
    using container_type::end;
    using container_type::cend;

    // Implement what is needed to be a Vector_i
    size_type size() const override { return container_type::size(); }
    size_type n_elem() const override { return size(); }
    Scalar operator()(size_type i) const override { return (*this)[i]; }
    Scalar operator[](size_type i) const override {
        return container_type::operator[](i);
    }

    // Some extra stuff we need for testing
    Scalar& operator()(size_type i) { return (*this)[i]; }
    Scalar& operator[](size_type i) { return container_type::operator[](i); }

    VectorModel(std::vector<Scalar> v) : container_type(v) {}
    VectorModel(size_type c, bool initialise) : container_type(c) {
        if (initialise) {
            for (auto& elem : *this) elem = 0;
        }
    }
};

template <typename Vector>
class TestingLibrary {
  public:
    typedef Vector vector_type;
    typedef typename vector_type::size_type size_type;
    typedef typename vector_type::scalar_type scalar_type;
    typedef typename vector_type::real_type real_type;

    /** Construct a testing library instance in order to check
     *  all basic functionality of a stored vector
     *
     *  \param prefix A prefix to use in the rapidcheck description string
     */
    TestingLibrary(std::string prefix = "")
          : m_prefix{prefix}, m_gen{argsgen} {}

    //, m_gen{argsgen, identity, identity} {}

    void run_checks() const;

  private:
    // The testing library and caller type
    typedef std::vector<scalar_type> data_type;
    typedef VectorModel<scalar_type> model_type;
    typedef vector_tests::ComparativeTests<model_type, vector_type> comptests;
    typedef RCTestableGenerator<model_type, vector_type, data_type> gen_type;

    static constexpr bool cplx = krims::IsComplexNumber<scalar_type>::value;
    static constexpr double genscale = cplx ? 0.8 : 1.0;

    /** Argument generation */
    static constexpr data_type argsgen() {
        return *gen::scale(genscale,
                           gen::numeric_container<std::vector<scalar_type>>())
                      .as("Vector data");
    }

    void once_test_initialiser_list_constructor() const;
    void once_test_as_scalar() const;

    std::string m_prefix;
    gen_type m_gen;
};

//
// ---------------------------------------
//

template <typename Matrix>
void TestingLibrary<Matrix>::once_test_initialiser_list_constructor() const {
    vector_type v{11.0, 22.0, 33.0, 44.0};

    CHECK((v.n_elem() == 4));
    CHECK((v.size() == 4));

    for (size_type i = 0; i < v.size(); ++i) {
        CHECK((v(i) == 11. * (i + 1)));
        CHECK((v[i] == 11. * (i + 1)));
    }

    vector_type vref(4);
    vref[0] = 11.0;
    vref[1] = 22.0;
    vref[2] = 33.0;
    vref[3] = 44.0;
    CHECK(v == vref);
}

template <typename Matrix>
void TestingLibrary<Matrix>::once_test_as_scalar() const {
    // TODO convert into standard test in indexable_tests.hh

    vector_type v{42.};
    scalar_type s = as_scalar(v);
    CHECK(s == v[0]);
}

template <typename Matrix>
void TestingLibrary<Matrix>::run_checks() const {
    // Shorter aliases:
    const NumCompAccuracyLevel deflvl = NumCompAccuracyLevel::Default;
    const NumCompAccuracyLevel low = NumCompAccuracyLevel::Lower;
    const NumCompAccuracyLevel sloppy = NumCompAccuracyLevel::Sloppy;
    const NumCompAccuracyLevel supersloppy = NumCompAccuracyLevel::SuperSloppy;
    const NumCompAccuracyLevel eps = NumCompAccuracyLevel::MachinePrecision;
    const NumCompAccuracyLevel eps10 =
          NumCompAccuracyLevel::TenMachinePrecision;

    once_test_initialiser_list_constructor();
    once_test_as_scalar();

    // Copying and ==
    /* XXX
    CHECK(m_gen.run_test(m_prefix + "Test ==", comptests::test_equivalence,
                         eps));
                         */
    CHECK(m_gen.run_test(m_prefix + "Test copying stored vectors",
                         comptests::test_copy, eps));

    // Read-only element access:
    CHECK(m_gen.run_test(m_prefix + "Element access via []",
                         comptests::test_element_access_vectorised, eps));
    CHECK(m_gen.run_test(m_prefix + "Element access via ()",
                         comptests::test_element_access, eps));
    CHECK(m_gen.run_test(m_prefix + "Read-only iterator.",
                         comptests::test_readonly_iterator, eps));

    // Read-write element access
    CHECK(m_gen.run_test(m_prefix + "Altering elements via ()",
                         comptests::test_setting_elements, eps));
    CHECK(m_gen.run_test(m_prefix + "Altering elements via []",
                         comptests::test_setting_elements_vectorised, eps));
    CHECK(m_gen.run_test(m_prefix + "Altering elements via iterator",
                         comptests::test_readwrite_iterator, eps));

    // Standard operations and norms
    CHECK(m_gen.run_test(m_prefix + "l1 norm", comptests::test_norm_l1, eps10));
    CHECK(m_gen.run_test(m_prefix + "l2 norm", comptests::test_norm_l2, eps10));
    CHECK(m_gen.run_test(m_prefix + "linf norm", comptests::test_norm_linf,
                         eps10));
    CHECK(m_gen.run_test(m_prefix + "accumulate function",
                         comptests::test_accumulate,
                         cplx ? supersloppy : sloppy));
    CHECK(m_gen.run_test(m_prefix + "dot and cdot function with vector_type",
                         comptests::template test_dot<vector_type>,
                         cplx ? supersloppy : low));
    CHECK(m_gen.run_test(
          m_prefix + "dot and cdot function with real model vector",
          comptests::template test_dot<VectorModel<real_type>>,
          cplx ? sloppy : low));
    CHECK(m_gen.run_test(m_prefix + "elementwise functions",
                         comptests::test_elementwise));

    if (!cplx) {
        CHECK(m_gen.run_test(m_prefix + "min and max function",
                             comptests::test_minmax, eps));
    }

    // Operations
    CHECK(m_gen.run_test(m_prefix + "Multiplication by scalar",
                         comptests::test_multiply_scalar));
    CHECK(m_gen.run_test(m_prefix + "Divide by scalar",
                         comptests::test_divide_scalar,
                         cplx ? sloppy : deflvl));
    CHECK(m_gen.run_test(m_prefix + "Add a vector",
                         comptests::template test_add<vector_type>));
    CHECK(m_gen.run_test(m_prefix + "Subtract a vector",
                         comptests::template test_subtract<vector_type>));
}

}  // namespace stored_matrix_tests
}  // namespace tests
}  // namespace linalgwrap
