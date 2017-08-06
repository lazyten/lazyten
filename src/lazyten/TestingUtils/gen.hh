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
#include "gen/Numeric.hh"
#include "gen/NumericAround.hh"
#include "gen/NumericConstants.hh"
#include "gen/NumericContainer.hh"
#include "gen/NumericNonZero.hh"
#include "gen/NumericSize.hh"
#include "gen/NumericTensor.hh"
#include "gen/WithProperties.hh"
#include "gen/with_l2_norm.hh"
#include "gen/with_l2_norm_in_range.hh"

namespace lazyten {
namespace gen {
// Import rapidcheck gen into this namespace
// such that all generators are available via
// lazyten::gen::
using namespace rc::gen;
}  // namespace gen
}  // namespace lazyten
