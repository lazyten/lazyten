// Exception thrown when an iterative solver failed
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
#include <krims/ExceptionSystem.hh>
#include <sstream>
#include <string>

namespace lazyten {

/** Base exception for all iterative solver exceptions */
class SolverException : public krims::ExceptionBase {
 public:
  SolverException() { m_extra = "A solver failed to converge."; }
};

//
// ------------
//

/**
 * Define a solver exception with no parameter.
 */
#define DefSolverException0(Exception0, outsequence)     \
  class Exception0 : public ::lazyten::SolverException { \
   public:                                               \
    Exception0() {                                       \
      std::ostringstream ss;                             \
      ss outsequence;                                    \
      ::krims::ExceptionBase::m_extra = ss.str();        \
    }                                                    \
  }

/**
 * Define a solver exception with one parameter.
 */
#define DefSolverException1(Exception1, type1, name1, outsequence) \
  class Exception1 : public ::lazyten::SolverException {           \
   public:                                                         \
    Exception1(const type1 a1) : name1(a1) {                       \
      std::ostringstream ss;                                       \
      ss outsequence;                                              \
      ::krims::ExceptionBase::m_extra = ss.str();                  \
    }                                                              \
    const type1 name1;                                             \
  }

/**
 * Define a solver exception with two parameters.
 */
#define DefSolverException2(Exception2, type1, name1, type2, name2, outsequence) \
  class Exception2 : public ::lazyten::SolverException {                         \
   public:                                                                       \
    Exception2(const type1 a1, const type2 a2) : name1(a1), name2(a2) {          \
      std::ostringstream ss;                                                     \
      ss outsequence;                                                            \
      ::krims::ExceptionBase::m_extra = ss.str();                                \
    }                                                                            \
    const type1 name1;                                                           \
    const type2 name2;                                                           \
  }

/** \brief Assert that a certain condition is satisfied in an iterative
 * solver. If not fail the solver state and throw an exception.
 *
 * The fail message will be set to the message of the exception.
 */
#define solver_assert(condition, state, exception) \
  {                                                \
    if (!(condition)) {                            \
      auto __exc__cept = exception;                \
      state.fail(__exc__cept.extra());             \
      this->on_failed(state);                      \
                                                   \
      assert_throw(condition, exception);          \
    }                                              \
  }

//
// ---------------
//

/** Thrown if maximum number of iterations reached */
DefSolverException1(ExcMaximumNumberOfIterationsReached, size_t, max_iter,
                    << "Reached maximum number of iterations (" << max_iter
                    << ") in the iterative solver");

/** Thrown if solver parameters are wrong or inconsistent */
DefSolverException1(ExcInvalidSolverParametersEncountered, std::string, details,
                    << "The solver could not run, because the set of solver "
                       "control parameters is not valid. Details: "
                    << details);

}  // namespace lazyten
