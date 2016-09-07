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
#include "OstreamState.hh"
#include <linalgwrap/Constants.hh>
#include <linalgwrap/Matrix_i.hh>
#include <ostream>

namespace linalgwrap {

// Forward declaration
template <typename Scalar>
class Matrix_i;

namespace io {

/** Class to print a matrix to an ostream such that
 *  the maximum number of sensible information can be
 *  noted from it.
 *
 *  The idea is to use this class when printing to cout
 *  or cerr or similar, i.e. output streams where the
 *  user directly sees what is written.
 */
class MatrixPrinter {
    // This is a really quick and dirty implementation.
    //
    // For some ideas on how to do this more cleverly look at
    // /usr/include/armadillo_bits/arma_ostream_bones.hpp
    // /usr/include/armadillo_bits/arma_ostream_meat.hpp

  public:
    /** Default constructor */
    MatrixPrinter();

    /** \name Properties */
    ///@{
    /** \brief Set the maximal screen width assumed in this class
     *
     * We try to keep the number of columns used for printing smaller than
     * this value.
     * */
    void max_screen_width(std::streamsize w);

    /** \brief Get the maximal screen width assumed in this class
     * */
    std::streamsize max_screen_width() const;

    /** \brief Set the character prepended to a known zero.
     *
     * Via the notion of sparsity patterns and iterators the library knows
     * of the locations of some zero values. These are maked in the printed
     * string via this character.
     * */
    void known_zero_character(char c);

    /** \brief Get the character prepended to a known zero.
     * */
    char known_zero_character() const;
    ///@}

    /** Print a generic matrix to the ostream */
    template <typename Scalar>
    void print(const Matrix_i<Scalar>&, std::ostream& out) const;

  private:
    /** Determine the width to use for each element and
     * set the formatting flags of the ostream
     *
     * \param mat  The matrix to scan
     * */
    template <typename Scalar>
    std::streamsize determine_width(const Matrix_i<Scalar>& mat,
                                    std::ostream& out) const;

    /** Print a single element to the ostream, given an actual width
     *
     * \param elem    The element to print
     * \param width   The width to use
     * \param out     The output stream to print to
     * */
    template <typename Scalar>
    void print_elem(const Scalar& elem, std::streamsize width,
                    std::ostream& out) const;

    /** Print a zero, either one which is or is not explicitly stored
     *
     * \param stored   Is the zero explicitly stored?
     * */
    void print_zero(bool stored, std::streamsize width,
                    std::ostream& out) const;

    /** Assumed maximal screen width */
    std::streamsize m_max_screen_width;

    /** The character to append to 0 in order to indicate that
     *  this zero is not explicitly stored */
    char m_known_zero_char;
};

//
// ------------------------------------------------
//

inline MatrixPrinter::MatrixPrinter()
      : m_max_screen_width(80), m_known_zero_char('*') {}

inline void MatrixPrinter::max_screen_width(std::streamsize w) {
    m_max_screen_width = w;
}

inline std::streamsize MatrixPrinter::max_screen_width() const {
    return m_max_screen_width;
}

inline void MatrixPrinter::known_zero_character(char c) {
    m_known_zero_char = c;
}

inline char MatrixPrinter::known_zero_character() const {
    return m_known_zero_char;
}

template <typename Scalar>
void MatrixPrinter::print(const Matrix_i<Scalar>& mat,
                          std::ostream& out) const {
    typedef typename Matrix_i<Scalar>::size_type size_type;

    // TODO Truncate the number of rows and columns which are
    // printed once it gets too many

    assert_dbg(out, krims::ExcIO());

    // Save ostream state
    // Upon destruction of this class the state is automatically
    // restored to what it was before.
    OstreamState state(out);

    // Determine width:
    std::streamsize width = determine_width(mat, out);

    // Here we use the fact that iterators are designed to skip known zero
    // elements. Hence if we traverse the matrix using an iterator and
    // by looping over all matrix elements, we know for sure that those
    // we do not reach via the iterator are not stored!
    auto it = std::begin(mat);
    for (size_type row = 0; row < mat.n_rows(); ++row) {
        for (size_type col = 0; col < mat.n_cols(); ++col) {
            if (it.col() == col && it.row() == row) {
                // Element can be reached via the iterator
                // => not a known zero
                print_elem(*it, width, out);
                ++it;
            } else {
                // Element cannot be reached via the iterator
                // => known zero.
                print_zero(false, width, out);
            }
        }
        out << std::endl;
    }
}

template <typename Scalar>
std::streamsize MatrixPrinter::determine_width(const Matrix_i<Scalar>&,
                                               std::ostream& out) const {
    using namespace std;
    // TODO have different cases here. See armadillo for ideas
    // TODO consider m_max_screen_width

    // Change stream flags:
    out.unsetf(ios::showbase);
    out.unsetf(ios::uppercase);
    out.unsetf(ios::showpos);
    out.setf(ios::scientific);
    out.setf(ios::right);
    out.unsetf(ios::fixed);
    out.fill(' ');
    out.precision(4);

    // Calculate and set width:
    // 1   Empty space
    // 3   Minus, number and .
    // 4   precision
    // 5   e+?? and space
    streamsize width = 1 + 3 + 4 + 4;

    return width;
}

template <typename Scalar>
void MatrixPrinter::print_elem(const Scalar& elem, std::streamsize width,
                               std::ostream& out) const {
    // TODO Do something special for NaN and Inf

    if (elem == Constants<Scalar>::zero) {
        print_zero(true, width, out);
        return;
    }

    // For some reason we need to set this each time ...
    out.width(width);
    out << elem;
}

inline void MatrixPrinter::print_zero(bool stored, std::streamsize width,
                                      std::ostream& out) const {
    // For some reason we need to set this each time ...
    if (!stored) {
        // Build string:
        std::stringstream s;
        s << m_known_zero_char;
        s << "0";

        out.width(width);
        out << s.str();
    } else {
        out.width(width);
        out << "0";
    }
}

}  // namespace io
}  // namespace linalgwrap
