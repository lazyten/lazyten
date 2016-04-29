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
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/Matrix_i.hh"
#include <ostream>
#include <string>

namespace linalgwrap {
namespace io {

/** \brief Marker interface to inherit all file type classes from
 *
 * All classes inheriting from this class should provide proper implementations
 * of the methods we define here and furthermore the following static
 * attributes:
 *
 * - static const std::vector<std::string> extensions
 *   A list of extensions, which is commonly used with this file type.
 *   the first element at index 0 should be the default extension to add when
 *   writing such a file.
 *
 * Specialisations for writing some kinds of matrices specially (like complex
 * or sparse matrices) can also be provided when implementing this class.
 * */
class FileType_i {
  public:
    /** Write a labelled matrix to the ostream under the format represented by
     * this class.
     *  Use the provided label string to indicate the matrix */
    template <typename Scalar>
    void write(std::ostream&, const std::string&,
               const Matrix_i<Scalar>&) const {
        assert_dbg(false, ExcNotImplemented());
    }

    /** Write a non-labelled matrix to the ostream under the format represented
     * by this class. */
    template <typename Scalar>
    void write(std::ostream&, const Matrix_i<Scalar>&) const {
        assert_dbg(false, ExcNotImplemented());
    }

    /** Write a comment string **/
    virtual void write_comment(std::ostream&, const std::string&) const = 0;
};

}  // namespace io
}  // namespace linalgwrap
