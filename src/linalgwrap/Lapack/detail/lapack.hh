//
// Copyright (C) 2017 by the linalgwrap authors
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
#ifdef LINALGWRAP_HAVE_LAPACK
#include "LapackPackedMatrix.hh"

namespace linalgwrap {
namespace detail {

/** Run the generalised Lapack eigensolver dspgv
 *
 * The Ap array is copied in and destroyed internally
 * The Bp array is copied in as well
 *
 * \param Ap  A array in packed format (as generated from make_packed_lower)
 * \param Bp  B array in packed format (as generated from make_packed_lower)
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 *
 * \returns  The Choleski factorisation of B in lower-triangle packed format.
 */
LapackPackedMatrix<double> run_dspgv(LapackPackedMatrix<double> Ap,
                                     LapackPackedMatrix<double> Bp,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info);

/** Run the Lapack eigensolver dspev
 *
 * The Ap array is copied in and destroyed internally
 *
 * \param Ap  A array in packed format (as generated from make_packed_lower)
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 */
void run_dspev(LapackPackedMatrix<double> Ap, std::vector<double>& evals,
               std::vector<double>& evecs, int& info);

// TODO What about the many million other interesting lapack solvers ...
//
// E.g. dsygst instead of dspgst

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
