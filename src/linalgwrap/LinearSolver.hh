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
#include "Armadillo/ArmadilloMatrix.hh"
#include "Base/Solvers.hh"
#include <krims/GenMap.hh>
#include <memory>

namespace linalgwrap {

#ifndef LINALGWRAP_HAVE_ARMADILLO
static_assert(false, "We need armadillo for this at the moment.");
#endif

template <typename LinearProblem>
class LinearSolverState final : public LinearSolverStateBase<LinearProblem> {
  // Forward constructors:
  using LinearSolverStateBase<LinearProblem>::LinearSolverStateBase;
};

struct LinearSolverKeys final : public LinearSolverBaseKeys {
  /** The solver method to use. Type: string */
  static const std::string method;
};

/** \brief Envelope linear solver that calls some
 *          inner solver depending on certain criteria
 *
 * TODO extend doc
 *
 * \tparam LinearProblem   linear problem type
 */
template <typename LinearProblem>
class LinearSolver final : public LinearSolverBase<LinearSolverState<LinearProblem>> {
 public:
  //@{
  /** Forwarded types */
  typedef LinearSolverState<LinearProblem> state_type;
  typedef LinearSolverBase<state_type> base_type;
  typedef typename base_type::size_type size_type;
  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::real_type real_type;
  typedef typename base_type::multivector_type multivector_type;
  typedef typename base_type::const_multivector_type const_multivector_type;
  typedef typename base_type::stored_matrix_type stored_matrix_type;
  typedef typename base_type::matrix_type matrix_type;
  typedef typename base_type::linproblem_type linproblem_type;
  //@}
  //

  /** \name Constructor */
  //@{
  /** Construct a linear solver with the default parameters */
  LinearSolver() {}

  /** Construct an linear solver setting the parameters from the map */
  LinearSolver(const krims::GenMap& map) : LinearSolver() { update_control_params(map); }
  //@}

  /** \name Iteration control */
  ///@{

  //! The method to use to actually solve the underlying linear problem.
  std::string method = "auto";

  /** Update control parameters from GenMap */
  void update_control_params(const krims::GenMap& map) {
    base_type::update_control_params(map);
    method = map.at(LinearSolverKeys::method, method);
  }

  //  /** Get the current settings of all internal control parameters and
  //   *  update the GenMap accordingly.
  //   */
  //  void get_control_params(krims::GenMap& map) const {
  //    base_type::get_control_params(map);
  //    map.update(LinearSolverKeys::method, method);
  //  }
  ///@}

  virtual void solve_state(state_type& state) const;

 private:
  // TODO Dummy for now ... we use implicitly armadillo for everything.
  typedef typename multivector_type::vector_type vector_type;
  void solve_hermitian_problem_armadillo(const matrix_type& A, vector_type& x,
                                         const vector_type& b) const;
};

//
// -----------------------------------------------------------
//

template <typename LinearProblem>
void LinearSolver<LinearProblem>::solve_hermitian_problem_armadillo(
      const matrix_type& A, vector_type& x, const vector_type& b) const {
  static_assert(std::is_same<stored_matrix_type, ArmadilloMatrix<double>>::value,
                "This hack only works for armadillo matrices");

  static_assert(std::is_same<scalar_type, double>::value,
                "This simple implementation only works for double");

  assert_dbg(A.is_symmetric(), ExcMatrixNotSymmetric());
  assert_dbg(props_contained_in(OperatorProperties::RealSymmetric, A.properties()),
             ExcMatrixNotSymmetric());

  // Zero out solution vector:
  x.set_zero();

  // Make arma matrices and vectors
  const arma::Mat<double>& m_arma =
        as_stored(A).data();  //.t() is skipped since m is symmetric
  arma::Col<double> b_arma(b.memptr(), b.size());
  arma::Col<double> x_arma(x.memptr(), x.size(), false);

  // Solve the linear system:
  bool result = arma::solve(x_arma, m_arma, b_arma);
  assert_throw(result, SolverException());
}

template <typename LinearProblem>
void LinearSolver<LinearProblem>::solve_state(state_type& state) const {
  assert_dbg(!state.is_failed(), krims::ExcInvalidState("Cannot solve a failed state"));
  const linproblem_type& problem = state.problem();

  // TODO only real symmetric problems implemented atm
  assert_dbg(
        props_contained_in(OperatorProperties::RealSymmetric, problem.A().properties()),
        krims::ExcNotImplemented());

  for (size_t i = 0; i < problem.n_systems(); ++i) {
    solve_hermitian_problem_armadillo(problem.A(), state.solution()[i], problem.rhs()[i]);
  }
}

}  // namespace linalgwrap
