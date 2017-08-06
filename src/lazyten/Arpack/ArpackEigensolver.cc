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

#include "ArpackEigensolver.hh"

#ifdef LAZYTEN_HAVE_ARPACK
namespace lazyten {

const std::string ArpackEigensolverKeys::max_iter = "max_iter";
const std::string ArpackEigensolverKeys::n_arnoldi_vectors = "n_arnoldi_vectors";
const std::string ArpackEigensolverKeys::mode = "mode";

}  // namespace lazyten
#endif  // LAZYTEN_HAVE_ARPACK
