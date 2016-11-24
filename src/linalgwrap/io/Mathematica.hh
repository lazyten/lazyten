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
#include "FileType_i.hh"
#include <krims/TypeUtils.hh>
#include <regex>
#include <vector>

namespace linalgwrap {
namespace io {

class Mathematica : public FileType_i {
 public:
  /** Construct a Mathematica file type. */
  Mathematica(int precision = 16);

  /** Construct a Mathematica file type, which drops elements less than
   *  a given threshold */
  Mathematica(double thresh, int precision = 16);

  /** \brief Write a matrix to a stream in mathematica format.
   *
   * The matrix ends up being written as a Mathematica list, which can be
   * used under the label provided. Note, that mathematica does not allow
   * certain characters in the label, so be careful.
   */
  template <typename Scalar>
  void write(std::ostream& out, const std::string& label,
             const Matrix_i<Scalar>& mat) const;

  /** \brief Write a matrix to a stream in mathematica format.
   *
   * The matrix ends up beeing written as a 2D Mathematica list.
   */
  template <typename Scalar>
  void write(std::ostream& out, const Matrix_i<Scalar>& mat) const;

  /** Write a labelled multivector to the ostream under the format represented
   * by this class.
   *  Use the provided label string to indicate the multivector
   */
  template <typename Scalar>
  void write(std::ostream& out, const MultiVector<Vector_i<Scalar>>& mv) const;

  template <typename Scalar>
  void write(std::ostream& out, const std::string& label, Scalar s) const;

  template <typename Scalar>
  void write(std::ostream& out, Scalar s) const;

  /** Write a comment **/
  void write_comment(std::ostream& out, const std::string& comment) const override {
    out << "(* " << comment << " *)" << std::endl;
  }

  /** Sanitise a label string, such that it satisfies the requirements of
   * the FileType */
  std::string normalise_label(const std::string& label) const override {
    return std::regex_replace(label, std::regex("[^a-zA-Z0-9]"), std::string(""));
  }

  /** Check whether a label string is a valid mathematica label */
  virtual bool is_valid_label(const std::string& label) const override {
    // Mathematica labels (variable names) can only consist of
    // letters or numbers
    return std::regex_match(label, std::regex("[a-zA-Z0-9]*"));
  }

  /** Extensions Mathematica files typically use */
  static const std::vector<std::string> extensions;

 private:
  template <typename Scalar>
  void write_elem(std::ostream& out, Scalar value) const;

  //! Threshold below which elements are considered zero.
  double m_thresh;

  //! Should we check for value being less then eps before writing it?
  bool m_check_for_thresh;

  //! Precision to use
  int m_precision;
};

//
// ---------------------------------------
//

template <typename Scalar>
void Mathematica::write(std::ostream& out, const std::string& label,
                        const Matrix_i<Scalar>& mat) const {
  assert_dbg(is_valid_label(label),
             ExcInvalidDataForFileType(" Mathematica labels need to consist "
                                       "of only letters and numbers, so " +
                                       label + " is not valid."));

  assert_throw(out, krims::ExcIO());
  out << label << " = ";
  write(out, mat);
}

template <typename Scalar>
void Mathematica::write(std::ostream& out,
                        const MultiVector<Vector_i<Scalar>>& mv) const {
  typedef typename Vector_i<Scalar>::size_type size_type;
  assert_throw(out, krims::ExcIO());
  if (mv.n_vectors() == 0 || mv.n_elem() == 0) {
    out << std::endl;
    return;
  }

  out << "{";
  for (size_type row = 0; row < mv.n_elem(); ++row) {
    if (row != 0) out << "," << std::endl;
    out << "{";
    for (size_type col = 0; col < mv.n_vectors(); ++col) {
      if (col != 0) out << ",";
      write_elem(out, (mv[col])[row]);
    }
    out << "}";
  }
  out << "};" << std::endl;
}

template <typename Scalar>
void Mathematica::write(std::ostream& out, const Matrix_i<Scalar>& mat) const {
  typedef typename Matrix_i<Scalar>::size_type size_type;
  assert_throw(out, krims::ExcIO());
  if (mat.n_rows() == 0 || mat.n_cols() == 0) {
    out << std::endl;
    return;
  }

  out << "{";
  for (size_type row = 0; row < mat.n_rows(); ++row) {
    if (row != 0) out << "," << std::endl;
    out << "{";
    for (size_type col = 0; col < mat.n_cols(); ++col) {
      if (col != 0) out << ",";
      write_elem(out, mat(row, col));
    }
    out << "}";
  }
  out << "};" << std::endl;
}

template <typename Scalar>
void Mathematica::write(std::ostream& out, const std::string& label, Scalar s) const {
  assert_dbg(is_valid_label(label),
             ExcInvalidDataForFileType(" Mathematica labels need to consist "
                                       "of only letters and numbers, so " +
                                       label + " is not valid."));

  assert_throw(out, krims::ExcIO());
  out << label << " = ";
  write(out, s);
}

template <typename Scalar>
void Mathematica::write(std::ostream& out, Scalar s) const {
  write_elem(out, s);
  out << ";" << std::endl;
}

template <typename Scalar>
void Mathematica::write_elem(std::ostream& out, Scalar value) const {
  assert_dbg(out, krims::ExcIO());
  if (m_check_for_thresh && std::abs(value) < m_thresh) {
    out << 0;
  } else {
    // Write to string stream and replace "e" by the "*^"
    // Mathematica uses for scientific exponents.
    std::stringstream ss;
    ss << std::setprecision(m_precision) << value;
    std::string res = ss.str();
    size_t pos = res.find("e");
    if (pos != std::string::npos) {
      res.replace(pos, 1, "*^");
    }
    out << res;
  }
}

}  // namespace io
}  // namespace linalgwrap
