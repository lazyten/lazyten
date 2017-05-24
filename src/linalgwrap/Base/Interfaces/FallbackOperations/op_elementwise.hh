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
#include "detail/ApplyElementwise.hh"
#include "linalgwrap/Base/Interfaces/Indexable_i.hh"
#include "linalgwrap/Base/Interfaces/Vector_i.hh"
#include "macro_defs.hh"
#include <krims/Functionals.hh>

namespace linalgwrap {

/** Return an indexable object which represents the elementwise absolute value
 * of an indexable object */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, krims::AbsFctr> abs(
      const Indexable& i) {
  linalgwrap_called_fallback();
  return detail::ApplyElementwise<Indexable, krims::AbsFctr>{i, krims::AbsFctr{}};
}

/** Return an indexable object which represents the elementwise complex
 * conjugate of an indexable object. */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, krims::ConjFctr> conj(
      const Indexable& i) {
  linalgwrap_called_fallback();
  return detail::ApplyElementwise<Indexable, krims::ConjFctr>{i, krims::ConjFctr{}};
}

/** Return an indexable object which represents the elementwise square root */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, krims::SqrtFctr> sqrt(
      const Indexable& i) {
  linalgwrap_called_fallback();
  return detail::ApplyElementwise<Indexable, krims::SqrtFctr>{i, krims::SqrtFctr{}};
}

/** Return an indexable object which represents the elementwise square */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, krims::SquareFctr> square(
      const Indexable& i) {
  linalgwrap_called_fallback();
  return detail::ApplyElementwise<Indexable, krims::SquareFctr>{i, krims::SquareFctr{}};
}

}  // namespace linalgwrap
