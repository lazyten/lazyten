//
// Copyright (C) 2016-17 by the linalgwrap authors
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
#include "linalgwrap/MultiVector.hh"
#include "linalgwrap/TypeUtils.hh"
#include <krims/TypeUtils/CheaplyCopyable_i.hh>
#include <memory>

namespace linalgwrap {
/** Container of a solution to an eigenproblem
 *
 * \tparam Evector  The type used for the eigenvectors.
 * \tparam Evalue   The type used for the eigenvalues.
 */
template <typename Evalue, typename Evector>
struct Eigensolution : public krims::CheaplyCopyable_i {
  /** \name Type definitions */
  ///@{
  /** The size type used in this class */
  typedef typename Evector::size_type size_type;

  /** The type used for the eigenvectors */
  typedef Evector evector_type;

  /** The type used for the eigenvalues */
  typedef Evalue evalue_type;
  ///@}

  /* \name Member attributes */
  ///@{
  /** Pointer to the eigenvectors */
  std::shared_ptr<MultiVector<evector_type>> evectors_ptr;

  /** Pointer to the eigenvalues */
  std::shared_ptr<std::vector<evalue_type>> evalues_ptr;
  ///@}

  /** \name Convenience access to attributes and properties*/
  ///@{
  /** The number of eigenpairs */
  size_type n_ep() const { return evalues_ptr->size(); }

  /**  Access to the eigenvectors
   *
   * The ordering of the eigenvectors agrees with the ordering of the
   * eigenvalues
   * */
  MultiVector<evector_type>& evectors() { return *evectors_ptr; }

  /** Access to the eigenvectors non-const version. */
  const MultiVector<evector_type>& evectors() const { return *evectors_ptr; }

  /** Access to the eigenvalues
   *
   * The eigenvectors are by convention ordered by arithmetic value (real
   * eigenvalues)
   * or magnitude (complex eigenvalues)
   * */
  std::vector<evalue_type>& evalues() { return *evalues_ptr; }

  /** Access to the eigenvalues non-const version  */
  const std::vector<evalue_type>& evalues() const { return *evalues_ptr; }
  ///@}

  /** \brief Construct by using the provided pointers to store the
   *         data. */
  Eigensolution(std::shared_ptr<MultiVector<evector_type>> evectors_ptr_,
                std::shared_ptr<std::vector<evalue_type>> evalues_ptr_)
        : evectors_ptr(std::move(evectors_ptr_)), evalues_ptr(std::move(evalues_ptr_)) {
    assert_size(evalues_ptr->size(), evectors_ptr->n_vectors());
  }

  /** Default constructor: Allocate empty structures */
  Eigensolution()
        : evectors_ptr(new MultiVector<evector_type>{}),
          evalues_ptr(new std::vector<evalue_type>{}) {}

  //@{
  /** Default move assignment and destructors */
  ~Eigensolution() = default;
  Eigensolution(Eigensolution&&) = default;
  Eigensolution& operator=(Eigensolution&&) = default;
  //@}

  /** Copy constructor: Perform a *deep* copy of the data of other */
  Eigensolution(const Eigensolution& other)
        : evectors_ptr{new MultiVector<evector_type>(other.evectors())},
          evalues_ptr(new std::vector<evalue_type>(other.evalues())) {}

  /** Copy assignment: Perform a *deep* copy of the data of other */
  Eigensolution& operator=(const Eigensolution& other) {
    evectors_ptr.reset(new MultiVector<evector_type>(other.evectors()));
    evalues_ptr.reset(new std::vector<evalue_type>(other.evalues()));
    return *this;
  }
};

namespace detail {
/** Determine eigenvalue type for an eigenproblem
 *
 * \tparam isHermitian   Is the eigenproblem Hermitian
 * \tparam MatrixScalar  The scalar type of the problem matrix
 */
template <bool isHermitian, typename MatrixScalar>
using Eigenvalue_t_inner =
      typename std::conditional<isHermitian ||
                                      krims::IsComplexNumber<MatrixScalar>::value,
                                MatrixScalar, std::complex<MatrixScalar>>::type;
// If hermitian or complex, keep type, else make it complex.
// TODO in theory real eigenvalues should work in this case
// as well but to keep the interface simple, complex eigenvalues
// are used in this case as well.

template <bool isHermitian, typename Matrix>
using Eigenvalue_t = Eigenvalue_t_inner<isHermitian, typename Matrix::scalar_type>;

/** Determine eigenvector type for an eigenproblem
 *
 * \tparam isHermitian   Is the eigenproblem Hermitian
 * \tparam Matrix        The problem matrix
 */
template <bool isHermitian, typename Matrix>
using Eigenvector_t = typename StoredTypeOf<Matrix>::type::type_family::template vector<
      Eigenvalue_t<isHermitian, Matrix>>;
// Use vector corresponding to the underlying stored matrix type's
// type family but of the scalar type of the eigenvalue.
}  // namespace detail

/** Using statement to determine the eigensolution type for an
 *  Eigenproblem's isHermitian flag and problem matrix type.
 *
 * The type depends on whether we have a Hermitian
 * eigenproblem and whether we have a real or complex scalar
 * type. The precise logic is as follows:
 *
 * a) real matrix => complex eigenvalues and complex eigenvectors
 *
 * b) complex matrix => complex eigenvalues and complex eigenvectors
 *
 * c) real hermitian matrix => real eigenvalues and real eigenvectors
 *
 * d) complex hermitian matrix => complex eigenvalues and complex
 *                                eigenvectors
 * */
template <bool isHermitian, typename Matrix>
using EigensolutionTypeFor = Eigensolution<detail::Eigenvalue_t<isHermitian, Matrix>,
                                           detail::Eigenvector_t<isHermitian, Matrix>>;

}  // namespace linalgwrap
