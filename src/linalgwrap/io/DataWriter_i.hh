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
#include <linalgwrap/Matrix_i.hh>
#include <linalgwrap/Subscribable.hh>

namespace linalgwrap {
namespace io {

// TODO we need something hierachical like
//      writer.iteration("1").write("ekin", 2.4);
//      also: printlevel!
//
//      The levels below should be of the same writer type
//      (Tree structure)
//
//      if the write operations are implemented out of place, we
//      save the templating in Scalar!

/** Generic interface class to write data to a stream or file under in some
 *  kind of data format. The details are specified by the implementing
 *  classes
 *
 * \tparam Scalar   The scalar type to use for all data.
 */
template <typename Scalar>
class DataWriter_i : public Subscribable {
  public:
    /** Write a labelled matrix to the ostream under the format represented by
     * this class.
     *  Use the provided label string to indicate the matrix
     *
     * \return Is the writer still in a good state?
     *  */
    virtual bool write(const std::string& label,
                       const Matrix_i<Scalar>& mat) = 0;

    /** Write a non-labelled matrix to the stream repersented by this class
     * under the format represented by this class
     *
     * \return Is the writer still in a good state?
     * */
    virtual bool write(const Matrix_i<Scalar>& mat) = 0;

    /** Write a comment string
     *
     * \return Is the writer still in a good state?
     * **/
    virtual bool write_comment(const std::string&) = 0;

    /** Write an empty line
     *
     * \return Is the writer still in a good state?
     * */
    virtual bool write_empty_line() = 0;

    /** Write a string verbatim as it is.
     *
     * \return Is the writer still in a good state?
     * */
    virtual bool write_verbatim(const std::string& s) = 0;

    /** Is the writer still in a good state? */
    virtual bool good() const = 0;

    operator bool() const { return good(); }
};

}  // namespace io
}  // namespace linalgwrap