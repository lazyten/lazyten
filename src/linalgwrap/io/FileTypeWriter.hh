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
#include "DataWriter_i.hh"
#include "FileType_i.hh"
#include <type_traits>

namespace linalgwrap {
namespace io {

/** \brief Write data to a std::ostream in a well-defined format.
 *
 * \tparam Scalar   The scalar type to use for all data.
 * \tparam FileType The file type to use for the formatting of the
 *                  data.
 */
template <typename FileType, typename Scalar>
class FileTypeWriter : public DataWriter_i<Scalar> {
  static_assert(std::is_base_of<FileType_i, FileType>::value,
                "FileType must be a child of FileType_i.");

 public:
  typedef FileType file_type;

  /** Construct a formatted stream writer from an output stream and a
   *  filetype object. */
  FileTypeWriter(std::ostream& out, const file_type ft);

  /** Write a labelled matrix to the ostream under the format represented by
   * this class.
   *  Use the provided label string to indicate the matrix */
  bool write(const std::string& label, const Matrix_i<Scalar>& mat) override;

  /** Write a non-labelled matrix to the stream represented by this class
   * under the format represented by this class */
  bool write(const Matrix_i<Scalar>& mat) override;

  /** Write a labelled multivector to the ostream under the format represented
   * by this class.
   *  Use the provided label string to indicate the multivector
   */
  bool write(const std::string& label,
             const MultiVector<Vector_i<Scalar>>& vecs) override;

  /** Write a non-labelled multivector to the stream represented by this class
   * under the format represented by this class */
  bool write(const MultiVector<Vector_i<Scalar>>& vecs) override;

  /** Write a labelled scalar under the format represented by this class
   * \return Is the writer still in a good state?
   * */
  virtual bool write(const std::string& label, Scalar s) override;

  /** Write a non-labelled scalar under the format represented by this class
   * \return Is the writer still in a good state?
   * */
  virtual bool write(Scalar s) override;

  /** Write a comment string **/
  bool write_comment(const std::string&) override;

  /** Write an empty line */
  bool write_empty_line() override;

  /** Write a string verbatim as it is. */
  bool write_verbatim(const std::string& s) override;

  /** Is the writer still in a good state? */
  bool good() const override;

  /** Sanitise a label string, such that it satisfies the required format */
  std::string normalise_label(const std::string& label) const override {
    return m_ft.normalise_label(label);
  }

  /** Check whether a label string is in the required format.*/
  bool is_valid_label(const std::string& label) const override {
    return m_ft.is_valid_label(label);
  }

 private:
  std::ostream& m_out;
  const file_type m_ft;
};

/** Convenience function to make a FormattedStreamWriter object by
 * constructing the FileType in-place */
template <typename FileType, typename Scalar, typename... GArgs>
FileTypeWriter<FileType, Scalar> make_writer(std::ostream& out, GArgs&&... gargs);

//
// ------------------------------------------------------------------------
//

template <typename FileType, typename Scalar>
inline FileTypeWriter<FileType, Scalar>::FileTypeWriter(std::ostream& out,
                                                        const file_type ft)
      : m_out{out}, m_ft{std::forward<const file_type>(ft)} {
  assert_throw(m_out, krims::ExcIO());
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::write(const std::string& label,
                                                    const Matrix_i<Scalar>& mat) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, label, mat);
  return m_out.good();
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::write(const Matrix_i<Scalar>& mat) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, mat);
  return m_out.good();
}

template <typename FileType, typename Scalar>
bool FileTypeWriter<FileType, Scalar>::write(const std::string& label,
                                             const MultiVector<Vector_i<Scalar>>& vecs) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, label, vecs);
  return m_out.good();
}

template <typename FileType, typename Scalar>
bool FileTypeWriter<FileType, Scalar>::write(const MultiVector<Vector_i<Scalar>>& vecs) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, vecs);
  return m_out.good();
}

template <typename FileType, typename Scalar>
bool FileTypeWriter<FileType, Scalar>::write(const std::string& label, Scalar s) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, label, s);
  return m_out.good();
}

template <typename FileType, typename Scalar>
bool FileTypeWriter<FileType, Scalar>::write(Scalar s) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write(m_out, s);
  return m_out.good();
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::write_comment(const std::string& comment) {
  assert_throw(m_out, krims::ExcIO());
  m_ft.write_comment(m_out, comment);
  return m_out.good();
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::write_empty_line() {
  assert_throw(m_out, krims::ExcIO());
  m_out << std::endl;
  return m_out.good();
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::write_verbatim(const std::string& s) {
  assert_throw(m_out, krims::ExcIO());
  m_out << s;
  return m_out.good();
}

template <typename FileType, typename Scalar>
inline bool FileTypeWriter<FileType, Scalar>::good() const {
  return m_out.good();
}

template <typename FileType, typename Scalar, typename... GArgs>
FileTypeWriter<FileType, Scalar> make_writer(std::ostream& out, GArgs&&... gargs) {
  return FileTypeWriter<FileType, Scalar>{out, FileType{gargs...}};
}

}  // namespace io
}  // namespace linalgwrap
