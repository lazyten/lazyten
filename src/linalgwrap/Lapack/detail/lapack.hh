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
#include "linalgwrap/config.hh"
#ifdef LINALGWRAP_HAVE_LAPACK

#include "LapackPackedMatrix.hh"
#include "LapackSymmetricMatrix.hh"

namespace linalgwrap {
namespace detail {

/** Run the Lapack eigensolver dspev
 *
 * The ap array is copied in and destroyed internally
 *
 * \param ap  A array in packed format
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 */
void run_dspev(LapackPackedMatrix<double> ap, std::vector<double>& evals,
               std::vector<double>& evecs, int& info);

/** Run the Lapack eigensolver dsyev
 *
 * The A array is copied in and destroyed internally
 *
 * \param A  A array in symmetric lapack matrix format
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 */
void run_dsyev(LapackSymmetricMatrix<double> a, std::vector<double>& evals,
               std::vector<double>& evecs, int& info);

/** Run the generalised Lapack eigensolver dspgv
 *
 * The ap array is copied in and destroyed internally
 * The Bp array is copied in as well
 *
 * \param ap  A array in packed format
 * \param Bp  B array in packed format
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 *
 * \returns  The Choleski factorisation of B in lower-triangle packed format.
 */
LapackPackedMatrix<double> run_dspgv(LapackPackedMatrix<double> ap,
                                     LapackPackedMatrix<double> bp,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info);

/** Run the generalised Lapack eigensolver dsygv
 *
 * The A array is copied in and destroyed internally
 * The B array is copied in as well
 *
 * \param A  A array in symmetric lapack matrix format
 * \param B  B array in symmetric lapack matrix format
 * \param evecs   The eigenvectors (in Fortran format, i.e. column-major)
 *                (will be resized by the function)
 * \param evals   The eigenvalues ordered by value
 *                (will be resized by the function)
 * \param info   The info parameter returned by Lapack
 *
 * \returns  The Choleski factorisation of B in lower-triangle packed format.
 */
LapackPackedMatrix<double> run_dsygv(LapackSymmetricMatrix<double> a,
                                     LapackSymmetricMatrix<double> b,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info);

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
