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
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/Matrix_i.hh"
#include <ostream>
#include <string>

namespace linalgwrap {
namespace io {

/** Exception to indicate that the provided data is not valid for this FileType,
 * e.g.
 * since it violates the specification of the file type */
DefException1(ExcInvalidDataForFileType, std::string,
              << "Invalid data was passed for the file type:" << arg1);

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
  virtual ~FileType_i() = default;
  FileType_i() = default;
  FileType_i(const FileType_i&) = default;
  FileType_i(FileType_i&&) = default;
  FileType_i& operator=(FileType_i&&) = default;
  FileType_i& operator=(const FileType_i&) = default;

  /** Write a labelled matrix to the ostream under the format represented by
   * this class.
   *  Use the provided label string to indicate the matrix */
  template <typename Scalar>
  void write(std::ostream&, const std::string&, const Matrix_i<Scalar>&) const {
    assert_implemented(false);
  }

  /** Write a non-labelled matrix to the ostream under the format represented
   * by this class. */
  template <typename Scalar>
  void write(std::ostream&, const Matrix_i<Scalar>&) const {
    assert_implemented(false);
  }

  /** Write a labelled multivector to the ostream under the format represented
   * by this class.
   *  Use the provided label string to indicate the multivector
   */
  template <typename Scalar>
  void write(std::ostream&, const std::string&,
             const MultiVector<Vector_i<Scalar>>&) const {
    assert_implemented(false);
  }

  /** Write a non-labelled multivector to the stream represented by this class
   * under the format represented by this class */
  template <typename Scalar>
  void write(std::ostream&, const MultiVector<Vector_i<Scalar>>&) const {
    assert_implemented(false);
  }

  /** Write a labelled scalar under the format represented by this class
   * \return Is the writer still in a good state?
   * */
  template <typename Scalar>
  void write(std::ostream&, const std::string&, Scalar) const {
    assert_implemented(false);
  }

  /** Write a non-labelled scalar under the format represented by this class
   * \return Is the writer still in a good state?
   * */
  template <typename Scalar>
  void write(std::ostream&, const std::string&) const {
    assert_implemented(false);
  }

  /** Write a comment string **/
  virtual void write_comment(std::ostream&, const std::string&) const = 0;

  /** Sanitise a label string, such that it satisfies the requirements of
   * the FileType */
  virtual std::string normalise_label(const std::string& label) const { return label; }

  /** Check weather a label string satisfies the requirements of the FileType
   */
  virtual bool is_valid_label(const std::string&) const { return true; }
};

}  // namespace io
}  // namespace linalgwrap
