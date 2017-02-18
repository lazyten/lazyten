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
#include "EigensystemSolver.hh"
#include <iterator>
#include <linalgwrap/io.hh>

namespace linalgwrap {

/** Namespace of functions to diagnose problems and call to gather information in case
 * things fail.
 *
 * They usually produce fairly verbose output to debug files and/or an output stream
 * and are by far not optimised for speed. Their syntax may also differ a little from
 * the usual methods to solve problems in linalgwrap.
 *
 * Most methods do nothing in release mode.
 * */
namespace rescue {

namespace detail {
/** Get the name of the file to write the Debug data to
 *  and print it to the outstream
 **/
std::string debugout_file();

/** Get the filetypewriter to write debug data. Write messages to out. */
template <typename Scalar>
io::FileTypeWriter<io::Mathematica, Scalar>& debugout() {
  typedef io::Mathematica FileType;
  static std::ofstream file(debugout_file());
  static io::FileTypeWriter<FileType, Scalar> writer(file, FileType());
  return writer;
}
}  // namespace detail

// TODO  have various methods here that gather statistics about matrices
//       and attempt to solve problems in a way to gather and display.
//       Ideas:
//           - condition number

/** \brief Attempt a rescue solve for an eigenproblem in DEBUG mode
 *
 * If we are in DEBUG mode it attempts to solve for all eigenvalues
 * up to numerical precision and writes them to a std::cerr.
 * In case of an error or exception it writes an error to std::cerr
 * and returns. In RELEASE the function does nothing.
 *
 * \note It takes the usual GenMap \t eigensystem or \t EigensystemSolver
 * would also take, but some of the parameters (like tolerance or method)
 * are not honoured, but overwritten internally to make the rescue
 * solve work.
 **/
template <typename Eigenproblem>
void failed_eigenproblem(Eigenproblem problem, krims::GenMap params = krims::GenMap()) {
#ifdef DEBUG
  typedef typename Eigenproblem::real_type real_type;
  typedef typename Eigenproblem::scalar_type scalar_type;

  const std::string ind = "   - ";
  std::cerr << "Eigenproblem rescue:" << std::endl;

  if (problem.dim() > 1000) {
    std::cerr << ind
              << "Did not run further heuristics:   problem.dim() ==" << problem.dim()
              << "." << std::endl;
    return;
  }

  // Dump
  //
  std::cerr << ind << "Dumping eigenproblem in '" << detail::debugout_file() << "'"
            << std::endl;
  detail::debugout<scalar_type>().write("A", problem.A());
  if (Eigenproblem::generalised) {
    detail::debugout<scalar_type>().write("B", problem.B());
  }
  detail::debugout<scalar_type>().write("Diag", problem.Diag());

  // Solve
  //
  try {
    std::cerr << ind << "Attempting rescue eigensolve";

    // Overwrite some user parameters:
    params.update(EigensystemSolverKeys::method, "auto");
    params.update(EigensystemSolverKeys::tolerance,
                  linalgwrap::Constants<real_type>::default_tolerance);

    // Solve for full spectrum:
    problem.n_ep(linalgwrap::Constants<size_t>::all);
    const auto sol =
          EigensystemSolver<Eigenproblem>{params}.solve(problem).eigensolution();

    std::cerr << "  ...  full eigenspectrum of problem:" << std::endl << std::endl;
    std::cerr << "       ";
    std::ostream_iterator<scalar_type> out_it(std::cerr, "  ");
    std::copy(std::begin(sol.evalues()), std::end(sol.evalues()), out_it);
    std::cerr << std::endl << std::endl;
  } catch (...) {
    std::cerr << "   ...  but the attempt to solve for the full spectrum "
                 "failed as well."
              << std::endl;
  }
#endif  // DEBUG
}

}  // namespace rescue
}  // namespace linalgwrap
