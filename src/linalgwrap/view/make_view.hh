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
#include "linalgwrap/view/BlockView.hh"
#include "linalgwrap/view/ScaleView.hh"
#include "linalgwrap/view/TransposeView.hh"

/** \file view/make_view.hh
 *  \brief Global convenience functions in the view namespace to make views
 */

namespace linalgwrap {
namespace view {
namespace detail {
// Forward declare some views:
template <typename Matrix>
class TransposeView;

template <typename Matrix>
class ScaleView;

template <typename Matrix>
class BlockView;
}

/** Convenience function to make a BlockView
 *
 *\param m      The Matrix to take a block of
 *\param rows   The range of row indices to consider
 *\param cols   The range of column indices to consider
 */
template <typename Matrix>
detail::BlockView<Matrix> block(Matrix& m,
                                Range<typename Matrix::size_type> rows,
                                Range<typename Matrix::size_type> cols);

/** Convenience function to make a BlockView which selects
 *  a range of columns (i.e. a column view)
 *
 *\param m      The Matrix to take some columns of
 *\param cols   The range of columns to consider
 */
template <typename Matrix>
detail::BlockView<Matrix> columns(Matrix& m,
                                  Range<typename Matrix::size_type> cols);

/** Convenience function to make a BlockView which selects
 *  a range of rows (i.e. a row view)
 *
 *\param m      The Matrix to take some columns of
 *\param rows   The range of rows to consider
 */
template <typename Matrix>
detail::BlockView<Matrix> rows(Matrix& m,
                               Range<typename Matrix::size_type> rows);

/** Convenience function to make a ScaleView
 *
 *\param m   The Matrix to scale
 *\param s   The scaling coefficient to apply to all matrix entries
 */
template <typename Matrix>
detail::ScaleView<Matrix> scale(Matrix& m, typename Matrix::scalar_type s) {
    return detail::ScaleView<Matrix>(m, s);
}

/** Convenience function to make a TransposeView
 *
 *\param m   The Matrix to transpose
 */
template <typename Matrix>
detail::TransposeView<Matrix> transpose(Matrix& m);

//
// ---------------------------------------------------------------
//

template <typename Matrix>
detail::BlockView<Matrix> block(Matrix& m,
                                Range<typename Matrix::size_type> rows,
                                Range<typename Matrix::size_type> cols) {
    return detail::BlockView<Matrix>{m, rows, cols};
}

template <typename Matrix>
detail::BlockView<Matrix> columns(Matrix& m,
                                  Range<typename Matrix::size_type> cols) {
    return detail::BlockView<Matrix>{m, range(m.n_rows()), cols};
}

template <typename Matrix>
detail::BlockView<Matrix> rows(Matrix& m,
                               Range<typename Matrix::size_type> rows) {
    return detail::BlockView<Matrix>{m, rows, range(m.n_cols())};
}

template <typename Matrix>
detail::TransposeView<Matrix> transpose(Matrix& m) {
    return detail::TransposeView<Matrix>(m);
}

}  // namespace view
}  // namespace linalgwrap
