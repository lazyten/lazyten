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
#include "lazyten/Base/Interfaces/Indexable_i.hh"
#include "macro_defs.hh"
#include "op_elementwise.hh"
#include <numeric>

namespace lazyten {

/** \brief Sums of elementwise products */
///@{
/** Compute the sum of the elementwise product of the two indexable objects
 *
 * No complex conjugation is underdone.
 **/
template <typename Indexable1, typename Indexable2>
typename std::common_type<ValidIndexableScalarT<Indexable1>,
                          ValidIndexableScalarT<Indexable2>>::type
dot(const Indexable1& A, const Indexable2& B);

/** Compute the sum of the elementwise product of two indexable objects.
 *
 * Complex conjugation is underdone in the first argument A.
 */
template <typename Indexable1, typename Indexable2>
auto cdot(const Indexable1& A, const Indexable2& B) -> decltype(dot(conj(A), B)) {
  lazyten_called_fallback();
  return dot(conj(A), B);
}
///@}

//
// -----------------------------------------
//

template <typename Indexable1, typename Indexable2>
typename std::common_type<ValidIndexableScalarT<Indexable1>,
                          ValidIndexableScalarT<Indexable2>>::type
dot(const Indexable1& A, const Indexable2& B) {
  assert_size(A.n_elem(), B.n_elem());
  lazyten_called_fallback();
  typedef typename std::common_type<typename Indexable1::scalar_type,
                                    typename Indexable2::scalar_type>::type ret_type;

  return std::inner_product(A.begin(), A.end(), B.begin(), Constants<ret_type>::zero);
}

}  // namespace lazyten
