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
/** \file Includes the most basic inferfaces along with the
 * default implementation for all important operations */

/** Signal the calling of fallback functions by printing to cerr
 *  (useful for debugging whether the specialisations have worked) */
//#define LINALGWRAP_SIGNAL_FALLBACK

// Basic interfaces
#include "BaseInterfaces/Indexable_i.hh"
#include "BaseInterfaces/StoredVector_i.hh"
#include "BaseInterfaces/Vector_i.hh"

// Default implementations for important Vector operations.
#include "BaseInterfaces/FallbackOperations.hh"
