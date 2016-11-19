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

#include "ArpackEigensolver.hh"

#ifdef LINALGWRAP_HAVE_ARPACK
namespace linalgwrap {

const std::string ArpackEigensolverKeys::max_iter = "max_iter";
const std::string ArpackEigensolverKeys::n_arnoldi_vectors =
      "n_arnoldi_vectors";
const std::string ArpackEigensolverKeys::mode = "mode";

}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_ARPACK
