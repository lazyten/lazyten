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
#include <krims/ExceptionSystem.hh>
#include <linalgwrap/Base/Interfaces/OperatorProperties.hh>

namespace linalgwrap {
//
// Matrix properties
//

/** A matrix that was expected to be Hermetian is not */
DefExceptionMsg(ExcMatrixNotHermitian,
                "Encountered a non-Hermetian matrix where a Hermetian matrix "
                "was expected");

/** A matrix that was expected to be Hermetian is not */
DefExceptionMsg(ExcMatrixNotSymmetric,
                "Encountered a non-Symmetric matrix where a Symmetric matrix "
                "was expected");

/** A matrix which was expected to be dense, contained sparsity */
DefExceptionMsg(ExcMatrixNotDense,
                "Encountered a non-dense matrix where a dense matrix was expected");

/** A matrix which was expected to be quadratic, was not */
DefExceptionMsg(ExcMatrixNotSquare,
                "Encountered a non-square matrix where a square matrix was expected");

/** Exception to indicate that an expected operator property is not satisfied
 *  by a matrix or operator.
 *
 *  Essentially a more generic case to the exceptions above. It should only be used
 *  in cases which are so general that throwing the exceptions above is difficult.
 */
DefException1(ExcOperatorPropertiesNotSatisfied, OperatorProperties,
              << "The matrix does not satisfy the added properties: " << arg1);

}  // linalgwrap
