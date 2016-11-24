// Exception thrown when an iterative solver failed
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
#include <krims/ExceptionSystem.hh>
#include <sstream>
#include <string>

namespace linalgwrap {

/** Base exception for all iterative solver exceptions */
class SolverException : public krims::ExceptionBase {
 public:
  /** Print exception-specific extra information to the outstream */
  void print_extra(std::ostream& out) const noexcept override {
    out << "A solver failed to converge.";
  }
};

//
// ------------
//

/**
 * Define a solver exception with no parameter.
 */
#define DefSolverException0(Exception0, outsequence)                      \
  class Exception0 : public ::linalgwrap::SolverException {               \
   public:                                                                \
    Exception0() {}                                                       \
    virtual ~Exception0() noexcept {}                                     \
    virtual void print_extra(std::ostream& out) const noexcept override { \
      out outsequence << std::endl;                                       \
    }                                                                     \
  }

/**
 * Define a solver exception with one parameter.
 */
#define DefSolverException1(Exception1, type1, name1, outsequence)        \
  class Exception1 : public ::linalgwrap::SolverException {               \
   public:                                                                \
    Exception1(const type1 a1) : name1(a1) {}                             \
    virtual ~Exception1() noexcept {}                                     \
    virtual void print_extra(std::ostream& out) const noexcept override { \
      out outsequence << std::endl;                                       \
    }                                                                     \
    const type1 name1;                                                    \
  }

/**
 * Define a solver exception with two parameters.
 */
#define DefSolverException2(Exception2, type1, name1, type2, name2, outsequence) \
  class Exception2 : public ::linalgwrap::SolverException {                      \
   public:                                                                       \
    Exception2(const type1 a1, const type2 a2) : name1(a1), name2(a2) {}         \
    virtual ~Exception2() noexcept {}                                            \
    virtual void print_extra(std::ostream& out) const noexcept override {        \
      out outsequence << std::endl;                                              \
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
      std::stringstream __sTriNg_str;              \
      auto __exc__cept = exception;                \
      __exc__cept.print_extra(__sTriNg_str);       \
                                                   \
      state.fail(__sTriNg_str.str());              \
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

}  // namespace linalgwrap
