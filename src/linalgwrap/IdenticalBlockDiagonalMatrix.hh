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
#include "BlockDiagonalMatrix.hh"

namespace linalgwrap {

/** Translation class to translate the IdenticalBlockDiagonalMatrix interface
* to one suitable to inherit off BlockDiagonalMatrixBase */
template <size_t N, typename... Matrices>
class IdenticalBlockDiagonalMatrixTParamRepeat;

//! Specialisation to terminate the buildup recursion
template <typename... Matrices>
class IdenticalBlockDiagonalMatrixTParamRepeat<0, Matrices...>
      : public BlockDiagonalMatrixBase<Matrices...> {};

//! Specialise to start the buildup recursion
template <size_t N, typename Matrix, typename... Matrices>
class IdenticalBlockDiagonalMatrixTParamRepeat<N, Matrix, Matrices...>
      : public IdenticalBlockDiagonalMatrixTParamRepeat<N - 1, Matrix, Matrix,
                                                        Matrices...> {};

/** Block diagonal matrix with all blocks identical */
template <size_t N, typename Matrix>
class IdenticalBlockDiagonalMatrix
      : public IdenticalBlockDiagonalMatrixTParamRepeat<N, Matrix> {};

// TODO check that block is quadratic!
}
