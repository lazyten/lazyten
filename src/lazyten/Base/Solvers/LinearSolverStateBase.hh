//
// Copyright (C) 2017 by the lazyten authors
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
#include "LinearProblem.hh"
#include "SolverStateBase.hh"

namespace lazyten {

/** \brief Class for the basic state of a linear solver.
 *
 * \tparam Matrix        The type of the operator A in the problem $Ax = b$.
 * \tparam Vector        The type of the vectors $x$,$b$ in the problem $Ax = b$.
 */
template <typename LinearProblem>
class LinearSolverStateBase : public SolverStateBase {
 public:
  /** \name Type definitions */
  ///@{
  /* The type of the linear problem */
  typedef LinearProblem linproblem_type;

  /** The matrix type used within the linear problem */
  typedef typename linproblem_type::stored_matrix_type stored_matrix_type;

  /** The size type used in the matrices*/
  typedef typename linproblem_type::size_type size_type;

  /** The scalar type used in the matrices*/
  typedef typename linproblem_type::scalar_type scalar_type;

  /** The real type used in the matrices */
  typedef typename linproblem_type::real_type real_type;

  /** The vector type used for the vectors x in
   *  $A x = b$ */
  typedef typename linproblem_type::multivector_type multivector_type;

  /** The vector type used for the vectors b in
   *  $A x = b$ */
  typedef typename linproblem_type::const_multivector_type const_multivector_type;

  /** The type of the matrix A in $A x = \lambda x$ */
  typedef typename linproblem_type::matrix_type matrix_type;
  ///@}

  /** Access to the linear problem */
  const linproblem_type& problem() const { return m_linproblem; }

  /** The current estimate for the solution of the linear problem $Ax = b$,
   *  i.e. the vectors $x$ (Const version).
   */
  const multivector_type& solution() const { return *m_solution_ptr; }

  /** The current estimate for the solution of the linear problem $Ax = b$,
   * i.e. the vectors $x$. */
  multivector_type& solution() { return *m_solution_ptr; }

  /** Constructor for the linear solver state.
   *
   * \param problem  The linear problem to solve.
   * \param solution The container for the solution vectors $x$
   * */
  LinearSolverStateBase(const linproblem_type problem, multivector_type& solution)
        : m_linproblem(std::move(problem)),
          m_solution_ptr("LinearSolverStateBase", solution) {
    assert_size(m_linproblem.dim(), solution.n_elem());
    assert_size(m_linproblem.n_systems(), solution.n_vectors());
  }

  /** Setup the guess of this state by copying the relevant data from another state. */
  void obtain_guess_from(const LinearSolverStateBase& other) {
    // Copy solution if the objects are not identical.
    if (&other.solution() != &solution()) {
      solution() = other.solution();
    }
  }

  /** Setup the guess of this state by copying the guess solution inside */
  void obtain_guess_from(const multivector_type& other_soln) {
    // Copy solution if the objects are not identical.
    if (&solution() != &other_soln) {
      solution() = other_soln;
    }
  }

 private:
  const linproblem_type m_linproblem;

  /** solution of the problem */
  krims::SubscriptionPointer<multivector_type> m_solution_ptr;
};

//@{
/** \brief struct representing a type (std::true_type, std::false_type) which
 *  indicates whether T is derived from LinearSolverStateBase
 **/
template <typename T, typename = void>
struct IsLinearSolverState : public std::false_type {};

template <typename T>
struct IsLinearSolverState<T, krims::VoidType<typename T::linproblem_type>>
      : public std::is_base_of<LinearSolverStateBase<typename T::linproblem_type>, T> {};
//@}

}  // namespace lazyten
