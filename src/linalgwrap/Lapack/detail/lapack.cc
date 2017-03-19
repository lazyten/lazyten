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

#ifdef LINALGWRAP_HAVE_LAPACK
#include "lapack.hh"

namespace linalgwrap {
namespace detail {

//
// dspgv: Real generalized symmetric-definite eigenproblem (A and B in packed storage)
//
extern "C" void dspgv_(int* itype, char* jobz, char* uplo, int* n, double* ap, double* bp,
                       double* w, double* z, int* ldz, double* work, int* info);

LapackPackedMatrix<double> run_dspgv(LapackPackedMatrix<double> Ap,
                                     LapackPackedMatrix<double> Bp,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info) {
  assert_size(Ap.n, Bp.n);
  assert_size(Ap.elements.size(), Bp.elements.size());
  assert_size(Ap.n * (Ap.n + 1) / 2, Ap.elements.size());
  evecs.resize(Ap.n * Ap.n);
  evals.resize(Ap.n);

  int n = static_cast<int>(Ap.n);
  std::vector<double> work(3 * Ap.n);  //< The work type
  int itype = 1;                       //< Jobtype to do (here: A*v = \lambda*B*v
  char jobz = 'V';                     //< Compute eigenvalues and eigenvectors
  char uplo = 'L';                     //< Ap and Bp contain the lower triangles
  dspgv_(&itype, &jobz, &uplo, &n, Ap.elements.data(), Bp.elements.data(), evals.data(),
         evecs.data(), &n, work.data(), &info);

  // After the dspgv run Bp actually contains the choleski factorisation of Bp, which
  // might be of use for later.
  return Bp;
}

//
// dspev: Real non-general eigenproblem (A in packed storage)
//
extern "C" void dspev_(char* jobz, char* uplo, int* n, double* ap, double* w, double* z,
                       int* ldz, double* work, int* info);

void run_dspev(LapackPackedMatrix<double> Ap, std::vector<double>& evals,
               std::vector<double>& evecs, int& info) {
  assert_size(Ap.n * (Ap.n + 1) / 2, Ap.elements.size());
  evecs.resize(Ap.n * Ap.n);
  evals.resize(Ap.n);

  int n = static_cast<int>(Ap.n);
  std::vector<double> work(3 * Ap.n);  //< The work type
  char jobz = 'V';                     //< Compute eigenvalues and eigenvectors
  char uplo = 'L';                     //< Ap and Bp contain the lower triangles
  dspev_(&jobz, &uplo, &n, Ap.elements.data(), evals.data(), evecs.data(), &n,
         work.data(), &info);
}

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
