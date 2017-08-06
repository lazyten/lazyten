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
/** \file Includes the core interfaces and datastructures used
 *  by the iterative and direct solvers provided by the
 *  external libraries */

// Base solver
#include "Solvers/SolverBase.hh"
#include "Solvers/SolverExceptions.hh"
#include "Solvers/SolverStateBase.hh"

// Iterative methods base types and structs
#include "Solvers/IterativeStateWrapper.hh"
#include "Solvers/IterativeWrapper.hh"
#include "Solvers/IterativeWrapperKeys.hh"

// Linear solver base types and structs
#include "Solvers/LinearProblem.hh"
#include "Solvers/LinearSolverBase.hh"
#include "Solvers/LinearSolverBaseKeys.hh"
#include "Solvers/LinearSolverStateBase.hh"

// Eigensolvers base types and structs
#include "Solvers/Eigenproblem.hh"
#include "Solvers/Eigensolution.hh"
#include "Solvers/EigensolverBase.hh"
#include "Solvers/EigensolverBaseKeys.hh"
#include "Solvers/EigensolverStateBase.hh"

// Helper functions and utilities
#include "Solvers/select_eigenvalues.hh"
