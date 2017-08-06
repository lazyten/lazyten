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
#include "lazyten/config.hh"
#ifdef LAZYTEN_HAVE_ARPACK

#include <array>
#include <cstring>
#include <krims/ExceptionSystem.hh>
#include <memory>
#include <vector>

namespace lazyten {

//
// Implicitly restarted Lanczos method (IRLM)
//
// Perform iterations:
extern "C" void dsaupd_(int* ido, char* bmat, const unsigned int* n, char* which,
                        const unsigned int* nev, const double* tol, double* resid,
                        int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
                        double* workd, double* workl, int* lworkl, int* info);

// Obtain eigenvectors of IRLM
extern "C" void dseupd_(int* rvec, char* howmany, int* select, double* d, double* z,
                        int* ldz, double* sigma, char* bmat, const unsigned int* n,
                        char* which, const unsigned int* nev, const double* tol,
                        double* resid, int* ncv, double* v, int* ldv, int* iparam,
                        int* ipntr, double* workd, double* workl, int* lworkl, int* info);

//
// Implicitly restarted Arnoldi method (IRAM)
//
// Perform iterations:
extern "C" void dnaupd_(int* ido, char* bmat, const unsigned int* n, char* which,
                        const unsigned int* nev, const double* tol, double* resid,
                        int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
                        double* workd, double* workl, int* lworkl, int* info);

// Obtain eigenvectors of IRAM
extern "C" void dneupd_(int* rvec, char* howmany, int* select, double* dr, double* di,
                        double* z, int* ldz, double* sigmar, double* sigmai,
                        double* workev, char* bmat, const unsigned int* n, char* which,
                        const unsigned int* nev, const double* tol, double* resid,
                        int* ncv, double* v, int* ldv, int* iparam, int* ipntr,
                        double* workd, double* workl, int* lworkl, int* info);

namespace detail {

class ds_upd_wrapper {
 public:
  typedef double scalar_type;

  /** Call dsaupd function
   *
   * \note The working array for reverse communication can be found in the
   * workd member.
   *
   * \param  ido     The ido parameter of dsaupd
   * \param  iparam  The iparam parameter of dsaupd
   * \param  ipntr   The ipntr parameter of dsaupd
   * \param  info    The info parameter of dsaupd
   */
  void arnoldi_step(int& ido, std::array<int, 11>& iparam, int& info) {
    dsaupd_(&ido, bmat.data(), &n, which.data(), &nev, &tol, resid_ptr->data(), &ncv,
            v_ptr->data(), &ldv, iparam.data(), ipntr.data(), workd_ptr->data(),
            workl_ptr->data(), &lworkl, &info);
  }

  /** Call dseupd function
   *
   * \param sigma   The shift parameter sigma
   * \param iparam  The iparam parameter as passed to arnoldi_step
   * \param evalues The eigenvalues result array
   *                after the call the size will be exactly the number
   *                of computed eigenpairs
   * \param evectors The eigenvectors result array
   *                 after the call it will contain the eigenvectors
   *                 as a column-major array of size n*nev
   *                 (if compute_eigenvectors == false the array
   *                 is not touched)
   * \param info    The info parameter as returned from dseupd
   * \param compute_eigenvectors  Should eigenvectors be computed at all
   */
  void eigenpairs(double sigma, std::array<int, 11>& iparam, std::vector<double>& evalues,
                  std::vector<double>& evectors, int& info,
                  bool compute_eigenvectors = true) {
    int rvec = 0;  // Compute no eigenvectors

    // Setup required objects:
    char howmany = 'A';  // All

    std::vector<int> select(ncv, 1);
    evalues.resize(nev);

    int ldz = n;
    if (compute_eigenvectors) {
      rvec = 1;
      evectors.resize(n * nev);
    }

    dseupd_(&rvec, &howmany, select.data(), evalues.data(), evectors.data(), &ldz, &sigma,
            bmat.data(), &n, which.data(), &nev, &tol, resid_ptr->data(), &ncv,
            v_ptr->data(), &ldv, iparam.data(), ipntr.data(), workd_ptr->data(),
            workl_ptr->data(), &lworkl, &info);
  }

  /** The ipntr parameter, which points into important locations
   *  inside workd and workl */
  std::array<int, 14> ipntr;

  /** The v array pointer */
  std::shared_ptr<std::vector<scalar_type>> v_ptr;

  /** The workl array pointer */
  std::shared_ptr<std::vector<scalar_type>> workl_ptr;

  /**
   * \param generalised_  General eigenproblem?
   * \param dim_      Dimensionality of the eigenproblem
   * \param which_    Which eigenpairs to compute
   * \param n_eigenpairs_  How many to compute
   * \param tolerance_          Tolerance for the eigensolver
   * \param n_arnoldi_vectors_  Number of arnoldi vectors
   * \param resid_ptr_      Shared pointer to residual vector data.
   *                   Has to be of length dim_.
   * \param compute_eigenvectors  Should the eigenvectors be computed
   */
  ds_upd_wrapper(bool generalised_, size_t dim_, std::string which_, size_t n_eigenpairs_,
                 scalar_type tolerance_, size_t n_arnoldi_vectors_,
                 std::shared_ptr<std::vector<scalar_type>> resid_ptr_,
                 std::shared_ptr<std::vector<scalar_type>> workd_ptr_) {
    if (generalised_) {
      bmat = {{'G', '\0'}};
    } else {
      bmat = {{'I', '\0'}};
    }

    // Dimensionality of the problem
    n = static_cast<int>(dim_);

    if (which_ == "SR") {
      // Arpack does not understand SR if
      // we deal with a real problem
      // -> change to SA
      which = {{'S', 'A', '\0'}};
    } else if (which_ == "LR") {
      which = {{'L', 'A', '\0'}};
    } else {
      // Copy the which characters
      std::strncpy(which.data(), which_.c_str(), 2);
      which[2] = '\0';
    }

    // Number of eigenpairs to compute
    nev = static_cast<int>(n_eigenpairs_);

    // Copy the tolerance
    tol = tolerance_;

    // resid as copied from pointer.
    assert_size(dim_, resid_ptr_->size());
    resid_ptr = resid_ptr_;

    // ncv is the number of Arnoldi vectors
    // ldv is the leading dimension of v (pretty much n)
    ncv = static_cast<int>(n_arnoldi_vectors_);
    ldv = n;

    // Resize v: should be n * ncv
    v_ptr.reset(new std::vector<scalar_type>(ldv * ncv, 0.0));

    // Ipntr: zero the array:
    std::fill(ipntr.begin(), ipntr.end(), 0.0);

    // workd: Workspace of length 3n
    assert_size(3 * dim_, workd_ptr_->size());
    workd_ptr = workd_ptr_;

    // size of workl needs to be ncv^2+8*ncv
    lworkl = ncv * ncv + 8 * ncv;
    workl_ptr.reset(new std::vector<scalar_type>(lworkl, 0.0));
  }

 private:
  // dnaupd
  std::array<char, 2> bmat;
  unsigned int n;
  std::array<char, 3> which;
  unsigned int nev;
  scalar_type tol;
  std::shared_ptr<std::vector<scalar_type>> resid_ptr;
  int ncv;
  int ldv;
  std::shared_ptr<std::vector<scalar_type>> workd_ptr;
  int lworkl;
};

}  // namespace detail

}  // namespace lazyten
#endif
