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

#ifndef LINALG_LAZY_MATRIX_I_HPP_
#define LINALG_LAZY_MATRIX_I_HPP_

#include "linalgwrap/LazyMatrixExpression.hh"

namespace linalgwrap {
/** \brief Interface of lazy matrices
 *
 * All classes which implement this interface perform lazy evaluation, e.g. the
 * product is not directly performed, but much rather it is stored as a
 * lazy_matrix_product expression. Copying the classes should not copy the
 * actual data, but only references to the data (i.e. the classes should contain
 * only small amounts of data that lives on the stack).
 *
 * Assumptions the library makes when using the matrices implementing this
 * interface:
 *    - Copying these type of objects is cheap.
 *    - Multiplication and addition of these objects is associative
 *
 * From the abstact class ``Matrix_i`` we expect abstract lazy matrices to
 * implement:
 *   - ``size_type n_rows() const``
 *   - ``size_type n_cols() const``
 *   - ``scalar_type operator() (size_type row, size_type col) const``
 *
 * From ``LazyMatrixExpression`` we expect the implementation of
 *   - ``lazy_matrix_expression_ptr_type clone() const``
 * see the appropriate classes for details what these methods are required to
 * do.
 *
 * It is recommended to overload the ``operator*`` between the lazy matrix
 * and ``stored_matrix_type`` as well (as defied below) if a more performant
 * multiplication than plain element-by-element is possible.
 *
 * \tparam StoredMatrix   The type of stored matrix to use
 */
template <typename StoredMatrix>
class LazyMatrix_i : public LazyMatrixExpression<StoredMatrix> {

  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    //
    // Partial implementation of the interface of a LazyMatrixExpression
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    virtual void update(const krims::ParameterMap&) override {
        assert_dbg(false, krims::ExcNotImplemented());
    }

    /** A default implementation of the matrix-matrix product.
     *
     * This loops over all occupied elements of this matrix and
     * performs the product. Overload this in case there is a
     * more performant implementation for your matrix available.
     */
    stored_matrix_type operator*(const stored_matrix_type& in) const override;
};

template <typename StoredMatrix>
typename LazyMatrix_i<StoredMatrix>::stored_matrix_type LazyMatrix_i<
      StoredMatrix>::operator*(const stored_matrix_type& in) const {
    assert_size(this->n_cols(), in.n_rows());

    // Allocate return storage
    stored_matrix_type res(this->n_rows(), in.n_cols(), true);

    // Perform the product
    for (size_type k = 0; k < in.n_cols(); ++k) {
        for (auto it = this->begin(); it != this->end(); ++it) {
            res(it.row(), k) += *it * in(it.col(), k);
        }
    }

    return res;
}

}  // namespace linalg
#endif  // LINALG_MATRIX_I_HPP_
