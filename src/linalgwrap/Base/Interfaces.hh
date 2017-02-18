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
/** \file Includes the most basic inferfaces for vectors and matrices
 * along with the default implementation for all important operations */

/** Signal the calling of fallback functions by printing to cerr
 *  (useful for debugging whether the specialisations have worked) */
//#define LINALGWRAP_SIGNAL_FALLBACK

// Basic interfaces
#include "Interfaces/Indexable_i.hh"
#include "Interfaces/Stored_i.hh"

// Vectors
#include "Interfaces/IsStoredVector.hh"
#include "Interfaces/MutableMemoryVector_i.hh"
#include "Interfaces/MutableVector_i.hh"
#include "Interfaces/Vector_i.hh"

// Matrices
#include "Interfaces/OperatorProperties.hh"
#include "Interfaces/Transposed.hh"

// Default implementations for important Vector/Matrix/Indexable operations.
#include "Interfaces/FallbackOperations.hh"
