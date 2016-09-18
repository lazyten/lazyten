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
#include "linalgwrap/BaseInterfaces/Indexable_i.hh"
#include "linalgwrap/BaseInterfaces/Vector_i.hh"
#include "linalgwrap/detail/GenericFunctionals.hh"
#include "macro_defs.hh"

namespace linalgwrap {

/** Return an indexable object which represents the elementwise absolute value
 * of an indexable object */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, detail::AbsFctr> abs(
      const Indexable& i) {
    linalgwrap_called_fallback();
    return detail::ApplyElementwise<Indexable, detail::AbsFctr>{
          i, detail::AbsFctr{}};
}

/** Return an indexable object which represents the elementwise complex
 * conjugate of an indexable object. */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, detail::ConjFctr> conj(
      const Indexable& i) {
    linalgwrap_called_fallback();
    return detail::ApplyElementwise<Indexable, detail::ConjFctr>{
          i, detail::ConjFctr{}};
}

/** Return an indexable object which represents the elementwise square root */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, detail::SqrtFctr> sqrt(
      const Indexable& i) {
    linalgwrap_called_fallback();
    return detail::ApplyElementwise<Indexable, detail::SqrtFctr>{
          i, detail::SqrtFctr{}};
}

/** Return an indexable object which represents the elementwise square */
template <typename Indexable>
detail::ApplyElementwise<ValidIndexableT<Indexable>, detail::SquareFctr> square(
      const Indexable& i) {
    linalgwrap_called_fallback();
    return detail::ApplyElementwise<Indexable, detail::SquareFctr>{
          i, detail::SquareFctr{}};
}

}  // namespace linalgwrap
