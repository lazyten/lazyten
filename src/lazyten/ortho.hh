//
// Copyright (C) 2017 by the lazyten authors
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
#include "Matrix_i.hh"
#include "MultiVector.hh"
#include <krims/GenMap.hh>
#include <krims/TypeUtils/EnableIfLibrary.hh>

namespace lazyten {

DefException1(ExcInvalidInnerProduct, std::string, << arg1);

/** Orthonormalise a set of vectors with respect to an orthogonalisation function
 *
 * The inner product function is expected to take two arguments of the
 * type MultiVector<Vector> and return a matrix of values which are the
 * inner products of each vector with each other.
 **/
template <typename Vector, typename InnerProduct,
          typename = krims::enable_if_t<!IsMatrix<krims::decay_t<InnerProduct>>::value>>
MultiVector<Vector> ortho(const MultiVector<Vector>& vs, InnerProduct&& prod,
                          const krims::GenMap& /*params*/ = krims::GenMap()) {
  typedef typename Vector::scalar_type scalar_type;
  using lazyten::as_multivector;
  // TODO For now this is plain and stupid Gram-Schmidt
  //      Better use svd or Householder reflections once they are there.

  if (vs.n_vectors() == 0) return vs;

  // For now this is just a (slow) modified Gram-Schmidt procedure.
  MultiVector<Vector> done;
  done.reserve(vs.n_vectors());

  // Deal with first vector:
  done.push_back(Vector(vs[0]));
  auto norm0 = as_scalar(prod(as_multivector(vs[0]), as_multivector(vs[0])));
  assert_dbg(norm0 >= 0, ExcInvalidInnerProduct("<v|v> needs to be non-negative."));
  done[0] /= std::sqrt(norm0);

  for (size_t i = 1; i < vs.n_vectors(); ++i) {
    // Compute inner products with all vectors we had so far:
    auto innprods = prod(as_multivector(vs[i]), done);
    assert_internal(innprods.n_rows() == 1);
    assert_internal(innprods.n_cols() == i);

    done.push_back(Vector(vs[i]));
    assert_internal(done.n_vectors() == i + 1);
    for (size_t j = 0; j < i; ++j) {
      auto j_norm = as_scalar(prod(as_multivector(done[j]), as_multivector(done[j])));
      assert_dbg(j_norm >= 0, ExcInvalidInnerProduct("<v|v> needs to be non-negative."));
      done[i] -= innprods[j] * done[j] / j_norm;
    }
    scalar_type norm = as_scalar(prod(as_multivector(done[i]), as_multivector(done[i])));
    assert_dbg(norm >= 0, ExcInvalidInnerProduct("<v|v> needs to be non-negative."));
    done[i] /= std::sqrt(norm);
  }
  return done;
}

/** Orthonormalise a set of vectors with respect to a metric matrix,
 *  i.e. the goal of the function is to make the set of vectors M-orthogonal. */
template <typename Vector, typename Matrix,
          typename = krims::enable_if_t<IsMatrix<krims::decay_t<Matrix>>::value>>
MultiVector<Vector> ortho(const MultiVector<Vector>& vs, const Matrix& m,
                          const krims::GenMap& params = krims::GenMap()) {
  return ortho(vs, [&m](const MultiVector<const Vector>& v1,
                        const MultiVector<const Vector>& v2) { return dot(v1, m * v2); },
               params);
}

/** Orthonormalise a set of vectors with respect to the identity matrix */
template <typename Vector>
MultiVector<Vector> ortho(const MultiVector<Vector>& vs,
                          const krims::GenMap& params = krims::GenMap()) {
  return ortho(vs, [](const MultiVector<const Vector>& v1,
                      const MultiVector<const Vector>& v2) { return dot(v1, v2); },
               params);
}

}  // namespace lazyten
