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
#include <regex>
#include <vector>

namespace linalgwrap {
namespace io {

class Mathematica : public FileType_i {
  public:
    /** Construct a Mathematica file type. */
    Mathematica();

    /** Construct a Mathematica file type, which drops elements less than
     *  a given threshold */
    Mathematica(double thresh);

    /** \brief Write a matrix to a stream in mathematica format.
     *
     * The matrix ends up being written as a Mathematica list, which can be
     * used under the label provided. Note, that mathematica does not allow
     * certain characters in the label, so be careful.
     */
    template <typename Matrix>
    void write(std::ostream& out, const std::string& label,
               const Matrix& mat) const;

    /** \brief Write a matrix to a stream in mathematica format.
     *
     * The matrix ends up beeing written as a 2D Mathematica list.
     */
    template <typename Scalar>
    void write(std::ostream& out, const Matrix_i<Scalar>& mat) const;

    /** Write a comment **/
    void write_comment(std::ostream&, const std::string&) const override;

    /** Sanitise a label string, such that it satisfies the requirements of
     * the FileType */
    std::string normalise_label(const std::string& label) const override;

    /** Check whether a label string is a valid mathematica label */
    virtual bool is_valid_label(const std::string& label) const override;

    /** Extensions Mathematica files typically use */
    static const std::vector<std::string> extensions;

  private:
    template <typename Scalar>
    void write_elem(std::ostream& out, Scalar value) const;

    //! Threshold below which elements are considered zero.
    double m_thresh;

    //! Should we check for value being less then eps before writing it?
    bool m_check_for_thresh;
};

//
// ---------------------------------------
//

template <typename Matrix>
void Mathematica::write(std::ostream& out, const std::string& label,
                        const Matrix& mat) const {
    assert_dbg(is_valid_label(label),
               ExcInvalidDataForFileType(" Mathematica labels need to consist "
                                         "of only letters and numbers, so " +
                                         label + " is not valid."));

    assert_throw(out, ExcIO());
    out << label << " = ";
    write(out, mat);
}

template <typename Scalar>
void Mathematica::write(std::ostream& out, const Matrix_i<Scalar>& mat) const {
    assert_throw(out, ExcIO());
    if (mat.n_rows() == 0 || mat.n_cols() == 0) return;

    out << "{";
    for (auto row : range(mat.n_rows())) {
        if (row != 0) out << "," << std::endl;
        out << "{";
        for (auto col : range(mat.n_cols())) {
            if (col != 0) out << ",";
            write_elem(out, mat(row, col));
        }
        out << "}";
    }
    out << "};" << std::endl;
}

template <typename Scalar>
void Mathematica::write_elem(std::ostream& out, Scalar value) const {
    assert_dbg(out, ExcIO());
    if (m_check_for_thresh && std::abs(value) < m_thresh) {
        out << 0;
    } else {
        out << value;
    }
}

void Mathematica::write_comment(std::ostream& out,
                                const std::string& comment) const {
    out << "(* " << comment << " *)" << std::endl;
}

std::string Mathematica::normalise_label(const std::string& label) const {
    return std::regex_replace(label, std::regex("[^a-zA-Z0-9]"), "");
}

bool Mathematica::is_valid_label(const std::string& label) const {
    // Mathematica labels (variable names) can only consist of
    // letters or numbers
    return std::regex_match(label, std::regex("[a-zA-Z0-9]*"));
}

}  // namespace io
}  // namespace linalgwrap
