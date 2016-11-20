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

#include "generators.hh"
#include "rapidcheck_utils.hh"
#include <catch.hpp>
#include <linalgwrap/MultiVector.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/TestingUtils.hh>

namespace linalgwrap {
namespace tests {
using namespace rc;

namespace multi_vector_tests {
/** Generate a vector of the given size with vectors of the provided type */
template <typename Vector>
std::vector<Vector> gen_vectors(size_t n_vecs) {
    auto n_elem = *gen::numeric_size<2>().as("Number of elements per vector");
    return *gen::container<std::vector<Vector>>(
          n_vecs, gen::numeric_tensor<Vector>(n_elem));
}

/** Generate a vector of vectors of the provided type */
template <typename Vector>
std::vector<Vector> gen_vectors() {
    return gen_vectors<Vector>(*gen::numeric_size<2>().as("Number of vectors"));
}

/** Take the result of the gen_vectors function and make
 * a multivector out of it, taking ownership of some of the vectors
 * or none of them (own_none = true)
 */
template <typename Vector>
MultiVector<Vector> gen_multivector(std::vector<Vector>& vs,
                                    const bool own_none = false) {
    MultiVector<Vector> ret;

    for (auto& v : vs) {
        bool transfer_own = *gen::arbitrary<bool>();
        if (own_none || !transfer_own) {
            ret.push_back(v);
        } else {
            Vector copy(v);
            ret.push_back(std::move(copy));
        }
    }

    return ret;
}
}  // namespace multi_vector_tests

TEST_CASE("MultiVector class", "[MultiVector]") {
    using namespace multi_vector_tests;

    typedef double scalar_type;
    typedef SmallVector<scalar_type> vector_type;
    typedef typename vector_type::size_type size_type;

    SECTION("Basic constructors") {
        vector_type v{1, 2, 3, 4};
        vector_type u{42, 43, 44};
        vector_type ucopy(u);

        MultiVector<vector_type> mv1{v};
        REQUIRE(mv1.n_vectors() == 1u);
        REQUIRE(mv1.n_elem() == 4u);
        REQUIRE(mv1.at_ptr(0).is_shared_ptr() == false);
        REQUIRE(mv1.at_ptr(0).get() == &v);
        REQUIRE(mv1[0] == v);

        MultiVector<vector_type> mv2{std::move(ucopy)};
        REQUIRE(mv2.n_vectors() == 1u);
        REQUIRE(mv2.n_elem() == 3u);
        REQUIRE(mv2.at_ptr(0).is_shared_ptr() == true);
        REQUIRE(mv2[0] == u);

        MultiVector<vector_type> mv3(3, 4, true);
        REQUIRE(mv3.n_vectors() == 4);
        REQUIRE(mv3.n_elem() == 3);
        for (size_type i = 0; i < 4; ++i) {
            REQUIRE(mv3.at_ptr(i).is_shared_ptr() == true);
            for (size_type j = 0; j < 3; ++j) {
                REQUIRE((mv3[i])[j] == 0);
            }
        }

        MultiVector<const vector_type> mvconst(v);
        REQUIRE(mvconst.n_vectors() == 1u);
        REQUIRE(mvconst.n_elem() == 4u);
        REQUIRE(mvconst.at_ptr(0).is_shared_ptr() == false);
        REQUIRE(mvconst[0] == v);

        MultiVector<vector_type> mv4{{1, 2},   // elem 0
                                     {3, 4},   // elem 1
                                     {5, 6}};  // elem 2
        REQUIRE(mv4[0] == vector_type({1, 3, 5}));
        REQUIRE(mv4[1] == vector_type({2, 4, 6}));

    }  // Basic constructors

    SECTION("as_multivector() function") {
        vector_type v{1, 2, 3, 4};
        vector_type u{42, 43, 44};
        vector_type ucopy(u);

        auto mv1 = as_multivector(v);
        REQUIRE(mv1.n_vectors() == 1u);
        REQUIRE(mv1.n_elem() == 4u);
        REQUIRE(mv1.at_ptr(0).is_shared_ptr() == false);
        REQUIRE(mv1.at_ptr(0).get() == &v);
        REQUIRE(mv1[0] == v);

        auto mv2 = as_multivector(std::move(ucopy));
        REQUIRE(mv2.n_vectors() == 1u);
        REQUIRE(mv2.n_elem() == 3u);
        REQUIRE(mv2.at_ptr(0).is_shared_ptr() == true);
        REQUIRE(mv2[0] == u);
    }  // As multivec

    SECTION("Checking Multivector push_back") {
        vector_type v{4, 8, 7, 3, 0, 2};
        MultiVector<vector_type> mv{vector_type(v)};
        vector_type vmatch{1, 2, 3, 4, 5, 6};
        vector_type vsmall{1, 2, 3, 4};

        mv.push_back(vmatch);
        REQUIRE(mv.n_vectors() == 2);
        REQUIRE(mv[0] == v);
        REQUIRE(mv[1] == vmatch);
        REQUIRE(mv.at_ptr(1).is_shared_ptr() == false);

        vector_type vmove{6, 5, 4, 3, 2, 1};
        mv.push_back(std::move(vmove));
        REQUIRE(mv.n_vectors() == 3);
        REQUIRE(mv[0] == v);
        REQUIRE(mv[1] == vmatch);
        REQUIRE((mv[2] == vector_type{6, 5, 4, 3, 2, 1}));
        REQUIRE(mv.at_ptr(2).is_shared_ptr() == true);

#ifdef DEBUG
        REQUIRE_THROWS_AS(mv.push_back(vsmall), krims::ExcSizeMismatch);
#endif

        mv.resize(0);
        REQUIRE(mv.n_vectors() == 0);
    }  // Checking Multivector push_back

    SECTION("Checking Multivector vector access") {
        auto test = []() {
            auto vecs = gen_vectors<vector_type>();
            auto mv = gen_multivector(vecs);

            RC_ASSERT(mv.n_vectors() == vecs.size());
            if (vecs.size() > 0) {
                RC_ASSERT(mv.n_elem() == vecs[0].size());
                RC_ASSERT(mv.front() == vecs.front());
                RC_ASSERT(mv.back() == vecs.back());
            }

            auto itvecs = std::begin(vecs);
            auto itmv = std::begin(mv);
            for (; itvecs != std::end(vecs); ++itvecs, ++itmv) {
                RC_ASSERT(*itvecs == *itmv);
            }

            size_type i = 0;
            itmv = std::begin(mv);
            for (; itmv != std::end(mv); ++itmv, ++i) {
                RC_ASSERT(*itmv == vecs[i]);
                RC_ASSERT(*itmv == mv[i]);
                RC_ASSERT(*itmv == mv.at(i));
                RC_ASSERT(*itmv == *mv.at_ptr(i));
            }

            // Modify some elements
            size_type modifyv = *gen::inRange<size_type>(0u, mv.n_vectors());
            size_type modifyi = *gen::inRange<size_type>(0u, mv.n_elem());
            scalar_type val = *gen::arbitrary<scalar_type>();
            {  // via []
                MultiVector<vector_type> copy = mv.copy_deep();
                copy[modifyv][modifyi] = val;

                i = 0;
                itmv = std::begin(copy);
                for (; itmv != std::end(copy); ++itmv, ++i) {
                    RC_ASSERT(*itmv == copy[i]);
                    RC_ASSERT(*itmv == copy.at(i));
                    RC_ASSERT(*itmv == *copy.at_ptr(i));
                    if (i != modifyv) {
                        RC_ASSERT(*itmv == vecs[i]);
                    } else {
                        RC_ASSERT((*itmv)[modifyi] == val);
                    }
                }
            }

            {  // via iterator
                MultiVector<vector_type> copy = mv.copy_deep();
                (*(std::begin(copy) + modifyv))[modifyi] = val;

                i = 0;
                itmv = std::begin(copy);
                for (; itmv != std::end(copy); ++itmv, ++i) {
                    RC_ASSERT(*itmv == copy[i]);
                    RC_ASSERT(*itmv == copy.at(i));
                    RC_ASSERT(*itmv == *copy.at_ptr(i));
                    if (i != modifyv) {
                        RC_ASSERT(*itmv == vecs[i]);
                    } else {
                        RC_ASSERT((*itmv)[modifyi] == val);
                    }
                }
            }

            RC_ASSERT(mv.empty() == (vecs.size() == 0));
            mv.clear();
            RC_ASSERT(mv.empty());
        };

        CHECK(rc::check("Checking Multivector vector access", test));
    }  // Multivector vector access

    SECTION("Subview, deep and shallow copy") {
        auto test = []() {
            auto vecs = gen_vectors<vector_type>();
            auto mv = gen_multivector(vecs);

            // Deep copy
            auto copy = mv.copy_deep();
            RC_ASSERT(copy.n_vectors() == mv.n_vectors());
            RC_ASSERT(copy.n_elem() == mv.n_elem());
            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(copy.at_ptr(i).is_shared_ptr() == true);
                RC_ASSERT(copy[i] == vecs[i]);
            }

            // Shallow copy
            MultiVector<vector_type> view(mv);
            RC_ASSERT(view.n_vectors() == mv.n_vectors());
            RC_ASSERT(view.n_elem() == mv.n_elem());
            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(view.at_ptr(i).is_shared_ptr() ==
                          mv.at_ptr(i).is_shared_ptr());
                RC_ASSERT(view[i] == vecs[i]);
            }

            // subviews
            auto range = *gen::range_within<size_type>(0, vecs.size());
            MultiVector<const vector_type> csub = mv.csubview(range);
            MultiVector<vector_type> sub = mv.subview(range);
            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(view.at_ptr(i).is_shared_ptr() ==
                          mv.at_ptr(i).is_shared_ptr());
                RC_ASSERT(view[i] == vecs[i]);
            }
        };

        CHECK(rc::check("Subview, deep and shallow copy", test));
    }  // Subview, deep and shallow copy

    SECTION("Vector type conversion") {
        auto test = []() {
            auto vecs = gen_vectors<vector_type>();
            auto mv = gen_multivector(vecs);

            MultiVector<const vector_type> cmv = mv;
            MultiVector<MutableMemoryVector_i<scalar_type>> mmv = mv;

            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(cmv[i] == mv[i]);
                RC_ASSERT(mmv[i] == mv[i]);
            }
        };

        CHECK(rc::check("Vector type conversion", test));
    }  // Vector type conversion

    SECTION("Obtaining memptrs") {
        auto test = []() {
            auto vecs = gen_vectors<vector_type>();
            auto mv = gen_multivector(vecs, true);

            std::vector<scalar_type*> ptrs = mv.memptrs();
            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(vecs[i].memptr() == ptrs[i]);
            }

            MultiVector<const vector_type> cmv = mv;
            std::vector<const scalar_type*> cptrs = cmv.memptrs();
            for (size_type i = 0; i < vecs.size(); ++i) {
                RC_ASSERT(vecs[i].memptr() == cptrs[i]);
            }

            static_assert(std::is_same<decltype(*(cmv.memptrs()[0])),
                                       const scalar_type&>::value,
                          "The type returned for constant memptr access is not "
                          "constant");
        };
        CHECK(rc::check("MultiVector: Obtaining memptrs", test));
    }  // Obtaining memptrs

    SECTION("outer_prod_sum() on multivectors") {
        auto test = [] {
            auto vecs1 = gen_vectors<vector_type>();
            auto mv1 = gen_multivector(vecs1);

            auto vecs2 = gen_vectors<vector_type>(vecs1.size());
            auto mv2 = gen_multivector(vecs2);

            // compute some outer products
            auto res = outer_prod_sum(mv1, mv2);
            auto res2 = outer_prod_sum(mv1, mv1);

            RC_ASSERT(res.n_rows() == mv1.n_elem());
            RC_ASSERT(res.n_cols() == mv2.n_elem());
            RC_ASSERT(res2.n_cols() == mv1.n_elem());
            RC_ASSERT(res2.n_rows() == mv1.n_elem());
            RC_ASSERT(res2.is_symmetric());

            for (size_type i = 0; i < res.n_rows(); ++i) {
                for (size_type j = 0; j < res.n_cols(); ++j) {
                    scalar_type sum = 0;
                    for (size_type k = 0; k < mv1.n_vectors(); ++k) {
                        sum += mv1[k](i) * mv2[k](j);
                    }  // k
                    RC_ASSERT_NC(krims::numcomp(sum) == res(i, j));
                }  // j
            }      // i
        };

        CHECK(rc::check("MultiVector: outer_prod_sum()", test));
    }  // outer_sum() on multivectors

    // TODO full stateful test

}  // MultiVector class
}  // namespace tests
}  // namespace linalgwrap
