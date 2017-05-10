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
// dspev: Real non-general eigenproblem (A in packed storage)
//
extern "C" void dspev_(char* jobz, char* uplo, int* n, double* ap, double* w, double* z,
                       int* ldz, double* work, int* info);

void run_dspev(LapackPackedMatrix<double> ap, std::vector<double>& evals,
               std::vector<double>& evecs, int& info) {
  assert_size(ap.n * (ap.n + 1) / 2, ap.elements.size());
  evecs.resize(ap.n * ap.n);
  evals.resize(ap.n);

  int n = static_cast<int>(ap.n);
  std::vector<double> work(3 * ap.n);  //< The work type
  char jobz = 'V';                     //< Compute eigenvalues and eigenvectors
  char uplo = 'L';                     //< ap contain the lower triangles
  dspev_(&jobz, &uplo, &n, ap.elements.data(), evals.data(), evecs.data(), &n,
         work.data(), &info);
}

//
// dsyev: Real non-general eigenproblem (A as full matrix)
//
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
                       double* work, int* lwork, int* info);

void run_dsyev(LapackSymmetricMatrix<double> a, std::vector<double>& evals,
               std::vector<double>& evecs, int& info) {
  assert_size(a.n * a.n, a.elements.size());
  evals.resize(a.n);

  int n = static_cast<int>(a.n);
  char jobz = 'V';  //< Compute eigenvalues and eigenvectors
  char uplo = 'L';  //< Use the lower triangles of A and B

  // Determine optimal work array size:
  double wkopt;
  int lwork = -1;
  dsyev_(&jobz, &uplo, &n, nullptr, &n, nullptr, &wkopt, &lwork, &info);
  if (info != 0) return;  // Error!

  // Allocate work array:
  // check that we don't get a wrongfully large size and allocate work array
  const auto wksize = static_cast<size_t>(wkopt);
  assert_internal(wksize > 0 && wksize <= std::max<size_t>(a.n * a.n, 10000));
  std::vector<double> work(wksize);
  lwork = static_cast<int>(wksize);

  // Call to solve ...
  dsyev_(&jobz, &uplo, &n, a.elements.data(), &n, evals.data(), work.data(), &lwork,
         &info);

  // Copy eigenvectors (which are returned inside the A-array)
  evecs = std::move(a.elements);
}

//
// dspgv: Real generalized symmetric eigenproblem (A and B in packed storage)
//
extern "C" void dspgv_(int* itype, char* jobz, char* uplo, int* n, double* ap, double* bp,
                       double* w, double* z, int* ldz, double* work, int* info);

LapackPackedMatrix<double> run_dspgv(LapackPackedMatrix<double> ap,
                                     LapackPackedMatrix<double> bp,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info) {
  assert_size(ap.n, bp.n);
  assert_size(ap.elements.size(), bp.elements.size());
  assert_size(ap.n * (ap.n + 1) / 2, ap.elements.size());
  evecs.resize(ap.n * ap.n);
  evals.resize(ap.n);

  int n = static_cast<int>(ap.n);
  std::vector<double> work(3 * ap.n);  //< The work type
  int itype = 1;                       //< Jobtype to do (here: A*v = \lambda*B*v
  char jobz = 'V';                     //< Compute eigenvalues and eigenvectors
  char uplo = 'L';                     //< ap and bp contain the lower triangles
  dspgv_(&itype, &jobz, &uplo, &n, ap.elements.data(), bp.elements.data(), evals.data(),
         evecs.data(), &n, work.data(), &info);

  // After the dspgv run bp actually contains the choleski factorisation of bp, which
  // might be of use for later.
  return bp;
}

//
// dsygv: Real generalised symmetric eigenproblem (A and B as full matrices)
//
extern "C" void dsygv_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda,
                       double* b, int* ldb, double* w, double* work, int* lwork,
                       int* info);

LapackPackedMatrix<double> run_dsygv(LapackSymmetricMatrix<double> a,
                                     LapackSymmetricMatrix<double> b,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info) {
  assert_size(a.n, b.n);
  assert_size(a.elements.size(), b.elements.size());
  assert_size(a.n * a.n, a.elements.size());
  evals.resize(a.n);

  int n = static_cast<int>(a.n);
  int itype = 1;    //< Jobtype to do (here: A*v = \lambda*B*v
  char jobz = 'V';  //< Compute eigenvalues and eigenvectors
  char uplo = 'L';  //< Use the lower triangles of A and B

  // Determine optimal work array size:
  double wkopt;
  int lwork = -1;
  dsygv_(&itype, &jobz, &uplo, &n, nullptr, &n, nullptr, &n, nullptr, &wkopt, &lwork,
         &info);
  if (info != 0) return LapackPackedMatrix<double>{};  // Error!

  // check that we don't get a wrongfully large size and allocate work array
  const auto wksize = static_cast<size_t>(wkopt);
  assert_internal(wksize > 0 && wksize <= std::max<size_t>(a.n * a.n, 10000));
  std::vector<double> work(wksize);
  lwork = static_cast<int>(wksize);

  // Call to solve ...
  dsygv_(&itype, &jobz, &uplo, &n, a.elements.data(), &n, b.elements.data(), &n,
         evals.data(), work.data(), &lwork, &info);

  // Copy eigenvectors (which are returned inside the A-array)
  evecs = std::move(a.elements);

  // Build the matrix to return with the choleski decomposition of the B matrix
  // (This is stored in the first n_chol elements of the B.elements array)
  LapackPackedMatrix<double> ret;
  ret.n = b.n;
  const size_t n_chol = b.n * (b.n + 1) / 2;
  ret.elements.resize(n_chol);
  std::move(std::begin(b.elements),
            std::begin(b.elements) + static_cast<ptrdiff_t>(n_chol),
            std::begin(ret.elements));
  return ret;
}

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
