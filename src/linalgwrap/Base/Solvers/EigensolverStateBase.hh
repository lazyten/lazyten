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
#include "Eigensolution.hh"
#include "SolverStateBase.hh"

namespace linalgwrap {
/** \brief Eigensolver state base class
 *
 * Provides the most basic state of each Eigenstate, e.g.
 * the eigenproblem and the current guess for the eigensolution
 * (if available at this stage).
 *
 * It is assumed that all deriving Eigensolver states are
 * constructible from only the eigenproblem as first argument.
 *
 * \tparam Eigenproblem  The type of the eigenproblem
 */
template <typename Eigenproblem>
class EigensolverStateBase : public SolverStateBase {
 public:
  /** \name Type definitions */
  ///@{
  /** The type of the eigenproblem */
  typedef Eigenproblem eproblem_type;

  /** The stored matrix type of the eigenproblem */
  typedef typename eproblem_type::stored_matrix_type stored_matrix_type;

  /** The size type of the eigenproblem */
  typedef typename eproblem_type::size_type size_type;

  /** The scalar type of the eigenproblem
   *  (but not necessarily of the solution) */
  typedef typename eproblem_type::scalar_type scalar_type;

  /** The real type of the eigenproblem */
  typedef typename eproblem_type::real_type real_type;

  /** The type of the eigensolution object as determined by the
   * EigensolutionTypeFor using statement. **/
  typedef EigensolutionTypeFor<eproblem_type::hermitian,
                               typename eproblem_type::matrix_diag_type>
        esoln_type;

  /** The eigenvalue type we use.
   * \note For details how this type is constructed see the
   *       EigensolutionTypeFor using statement (in Eigensolution.hh).
   **/
  typedef typename esoln_type::evalue_type evalue_type;

  /** The eigenvector type we use
   * \note For details how this type is constructed see the Eigensolution_t
   *using statement.
   **/
  typedef typename esoln_type::evector_type evector_type;
  ///@}

  /** Access to eigenproblem */
  const eproblem_type& eigenproblem() const { return m_eigenproblem; }

  /** Const access to eigensolution */
  const esoln_type& eigensolution() const { return m_eigensolution; }

  /** Access to eigensolution*/
  esoln_type& eigensolution() { return m_eigensolution; }

  /** \brief Constructor
   *
   * Initialise the eigensolver state from an eigenproblem object,
   * which is copied into the current object
   * (This is fine since the eigenproblem object contains only
   * pointers or refererences)
   */
  EigensolverStateBase(const eproblem_type eigenproblem)
        : m_eigenproblem(std::move(eigenproblem)), m_eigensolution{} {};

 private:
  /* The eigenproblem we wish to solve (contains only pointers
   *  ->ok to store) */
  const eproblem_type m_eigenproblem;

  /* The container for the eigensolution (contains only pointers
   *  ->ok to store) */
  esoln_type m_eigensolution;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from EigensolverStateBase
 **/
template <typename T, typename = void>
struct IsEigensolverState : public std::false_type {};

template <typename T>
struct IsEigensolverState<T, krims::VoidType<typename T::eproblem_type>>
      : public std::is_base_of<EigensolverStateBase<typename T::eproblem_type>, T> {};
//@}

}  // namespace linalgwrap
