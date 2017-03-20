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

//
// dsyev: Real non-general eigenproblem (A as full matrix)
//
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w,
                       double* work, int* lwork, int* info);

void run_dsyev(LapackSymmetricMatrix<double> A, std::vector<double>& evals,
               std::vector<double>& evecs, int& info) {
  assert_size(A.n * A.n, A.elements.size());
  evals.resize(A.n);

  int n = static_cast<int>(A.n);
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
  assert_dbg(wksize > 0 && wksize <= std::max<size_t>(A.n * A.n, 1000),
             krims::ExcInternalError());
  std::vector<double> work(wksize);
  lwork = static_cast<int>(wksize);

  // Call to solve ...
  dsyev_(&jobz, &uplo, &n, A.elements.data(), &n, evals.data(), work.data(), &lwork,
         &info);

  // Copy eigenvectors (which are returned inside the A-array)
  evecs = std::move(A.elements);
}

//
// dspgv: Real generalized symmetric eigenproblem (A and B in packed storage)
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
// dsygv: Real generalised symmetric eigenproblem (A and B as full matrices)
//
extern "C" void dsygv_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda,
                       double* b, int* ldb, double* w, double* work, int* lwork,
                       int* info);

LapackPackedMatrix<double> run_dsygv(LapackSymmetricMatrix<double> A,
                                     LapackSymmetricMatrix<double> B,
                                     std::vector<double>& evals,
                                     std::vector<double>& evecs, int& info) {
  assert_size(A.n, B.n);
  assert_size(A.elements.size(), B.elements.size());
  assert_size(A.n * A.n, A.elements.size());
  evals.resize(A.n);

  int n = static_cast<int>(A.n);
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
  assert_dbg(wksize > 0 && wksize <= std::max<size_t>(A.n * A.n, 1000),
             krims::ExcInternalError());
  std::vector<double> work(wksize);
  lwork = static_cast<int>(wksize);

  // Call to solve ...
  dsygv_(&itype, &jobz, &uplo, &n, A.elements.data(), &n, B.elements.data(), &n,
         evals.data(), work.data(), &lwork, &info);

  // Copy eigenvectors (which are returned inside the A-array)
  evecs = std::move(A.elements);

  // Build the matrix to return with the choleski decomposition of the B matrix
  // (This is stored in the first n_chol elements of the B.elements array)
  LapackPackedMatrix<double> ret;
  ret.n = B.n;
  const size_t n_chol = B.n * (B.n + 1) / 2;
  ret.elements.resize(n_chol);
  std::move(std::begin(B.elements),
            std::begin(B.elements) + static_cast<ptrdiff_t>(n_chol),
            std::begin(ret.elements));
  return ret;
}

}  // namespace detail
}  // namespace linalgwrap
#endif  // LINALGWRAP_HAVE_LAPACK
