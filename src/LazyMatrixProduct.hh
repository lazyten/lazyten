#ifndef LINALG_MATRIX_PRODUCT_HPP_
#define LINALG_MATRIX_PRODUCT_HPP_

#include "ParameterMap.hh"
#include "LazyMatrixExpression.hh"
#include "Constants.hh"
#include <memory>
#include <vector>
#include <algorithm>
#include <iterator>

namespace linalgwrap {

// Forward declaration
template <typename StoredMatrix>
class LazyMatrixExpression;

/**
 * \brief Class to represent a product of matrices, which are scaled
 * by a common coefficient.
 */
template <typename StoredMatrix>
class LazyMatrixProduct : public LazyMatrixExpression<StoredMatrix> {
  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    /** A swap function for LazyMatrixProducts */
    friend void swap(LazyMatrixProduct& first, LazyMatrixProduct& second) {
        using std::swap;  // enable ADL

        swap(first.m_coefficient, second.m_coefficient);
        swap(first.m_factors, second.m_factors);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

    /** \name Constructors, desctructors and assignment
     */
    ///@{
    /** \brief Create a matrix product object:
     *
     * @param expr   The first matrix expression factor
     * @param factor The scalar factor to initialise the object with
     */
    explicit LazyMatrixProduct(const LazyMatrixExpression<StoredMatrix>& expr,
                               scalar_type factor = Constants<scalar_type>::one)
          : m_coefficient(factor) {

        // Push back the first factor:
        m_factors.push_back(std::move(expr.clone()));
    }

    /** \brief Copy constructor */
    LazyMatrixProduct(const LazyMatrixProduct&) = default;

    /** \brief Copy and scale constructor constructor
     *
     * Copy the contents of another product and scale it altogether
     * by an extra factor.
     *
     *  @param factor   The factor to scale the product with
     * */
    explicit LazyMatrixProduct(LazyMatrixProduct prod, scalar_type factor) {
        swap(*this, prod);
        m_coefficient *= factor;
    }

    /** \brief Default move constructor */
    LazyMatrixProduct(LazyMatrixProduct&&) = default;

    /** \brief Default destructor */
    ~LazyMatrixProduct() = default;

    /** Assignment operator */
    LazyMatrixProduct& operator=(LazyMatrixProduct other) {
        swap(*this, other);
        return *this;
    }
    ///@}

    //
    // Push back further factors
    //
    /** \brief Add a factor by copying or moving the matrix \p m into
     * the current object.
     */
    void push_factor(const LazyMatrixExpression<StoredMatrix>& e) {
        // Check fitting dimenisionality of matrix expressions:
        assert_size(n_cols(), e.n_rows());

        // Place into m_factors by moving a copy there
        m_factors.push_back(std::move(e.clone()));
    }

    /** \brief Push back all factors of a product onto another product
     */
    void push_factor(LazyMatrixProduct prod) {
        // Check fitting dimenisionality of matrix expressions:
        assert_size(n_cols(), prod.n_rows());

        // Move all factors of m to the very end:
        std::move(std::begin(prod.m_factors), std::end(prod.m_factors),
                  back_inserter(m_factors));

        // Adjust the scaling:
        m_coefficient *= prod.m_coefficient;
    }

    //
    // Other modifications of the product
    //
    /** \brief Scale the whole matrix_factor by the scalar \p c
     */
    void scale(const scalar_type c) {
        assert_finite(c);
        m_coefficient *= c;
    }

    //
    // Matrix_i interface
    //
    /** \brief Number of rows of the matrix */
    size_type n_rows() const override {
        // The number of rows of the product
        // equals the number of rows of the first element in the product.
        return m_factors.front()->n_rows();
    }

    /** \brief Number of columns of the matrix  */
    size_type n_cols() const override {
        // The number of columns of the product
        // equals the nuber of columns of the last element in the product.
        return m_factors.back()->n_cols();
    }

    /** \brief return an element of the matrix    */
    scalar_type operator()(size_type row, size_type col) const override;

    //
    // LazyMatrixExpression interface
    //
    /** \brief Extract a block of values out of the matrix and
     *         return it as a stored matrix of the appropriate size
     *
     * For more details of the interface see the function of the same
     * name in ``LazyMatrixExpression``.
     *
     * \param row_range   The Range object representing the range of rows
     *                    to extract. Note that it is a half-open interval
     *                    i.e. the LHS is inclusive, but the RHS not.
     * \param col_range   The Range object representing the range of
     *                    columns to extract.
     */
    virtual stored_matrix_type extract_block(
          Range<size_type> row_range,
          Range<size_type> col_range) const override;

    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap& map) override {
        // Pass the call onto all factors:
        for (auto& expression : m_factors) {
            expression->update(map);
        }
    }

    /** \brief Multiplication with a stored matrix */
    virtual stored_matrix_type operator*(
          const stored_matrix_type& m) const override {
        assert_size(n_cols(), m.n_rows());

        // Deal with last factor:
        auto it = m_factors.rbegin();
        stored_matrix_type res = m_coefficient * (*(*it) * m);
        ++it;

        // Deal with all others:
        contract_in_place(it, m_factors.rend(), res);
        return res;
    }

    /** \brief Print the expression tree to this outstream
     * */
    virtual void print_tree(std::ostream& o) const override {
        // TODO to be implemented
        assert_dbg(false, ExcNotImplemented());
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixProduct(*this));
    }

    //
    // In-place scalar operators:
    //
    /** \brief scale a matrix using a scalar
     */
    LazyMatrixProduct& operator*=(scalar_type c) {
        assert_finite(c);
        this->scale(c);
        return *this;
    }

    LazyMatrixProduct& operator/=(scalar_type c) {
        assert_finite(c);
        assert_dbg(c != 0, ExcDevideByZero());
        this->scale(Constants<scalar_type>::one / c);
        return *this;
    }

  private:
    /** \brief Contract a range of factors from m_factor with the stored_matrix
     *m
     * After the execution of the function \p m contains the result.
     *
     * The function is equivalent to
     * ```
     * for (; begin != end; ++begin) {
     *     m = (**begin) * m;
     * }
     */
    template <typename Iterator>
    void contract_in_place(Iterator begin, Iterator end,
                           stored_matrix_type& m) const;

    //! The vector containing the lazy matrix expression factors
    std::vector<lazy_matrix_expression_ptr_type> m_factors;
    scalar_type m_coefficient;
};

/** \brief Multiply two lazy matrix products */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(LazyMatrixProduct<StoredMatrix> lhs,
                                          LazyMatrixProduct<StoredMatrix> rhs) {
    lhs.push_factor(std::move(rhs));
    return lhs;
}

/** \brief Multiply a lazy matrix product and an expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      LazyMatrixProduct<StoredMatrix> lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {
    lhs.push_factor(rhs);
    return lhs;
}

/** \brief Multiply a lazy matrix product and an expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      const LazyMatrixExpression<StoredMatrix>& rhs,
      LazyMatrixProduct<StoredMatrix> lhs) {
    return lhs * rhs;
}

//
// Operator*
//
/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      LazyMatrixProduct<StoredMatrix> m, typename StoredMatrix::scalar_type s) {
    m *= s;
    return m;
}

/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(typename StoredMatrix::scalar_type s,
                                          LazyMatrixProduct<StoredMatrix> m) {
    return m * s;
}

/** \brief Scale a lazy matrix product  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator/(
      LazyMatrixProduct<StoredMatrix> lhs,
      typename StoredMatrix::scalar_type rhs) {
    lhs /= rhs;
    return lhs;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator-(LazyMatrixProduct<StoredMatrix> mat) {
    typedef typename StoredMatrix::scalar_type scalar_type;
    return -Constants<scalar_type>::one * mat;
}

//
// ------------------------------------------------------------
//
//
// LazyMatrixProduct
//
template <typename StoredMatrix>
typename LazyMatrixProduct<StoredMatrix>::scalar_type
      LazyMatrixProduct<StoredMatrix>::
      operator()(size_type row, size_type col) const {
    assert_range(0, row, n_rows());
    assert_range(0, col, n_cols());
    auto block = extract_block({row, row + 1}, {col, col + 1});
    return block(0, 0);
}

template <typename StoredMatrix>
typename LazyMatrixProduct<StoredMatrix>::stored_matrix_type LazyMatrixProduct<
      StoredMatrix>::extract_block(Range<size_type> row_range,
                                   Range<size_type> col_range) const {
    // At least one range is empty -> no work to be done:
    if (row_range.is_empty() || col_range.is_empty()) {
        return SmallMatrix<scalar_type>{row_range.length(), col_range.length()};
    }

    // Assertive checks:
    assert_upper_bound(row_range.last(), this->n_rows() + 1);
    assert_upper_bound(col_range.last(), this->n_cols() + 1);

    // If there is only one factor, just perform the operation downstream
    // and scale it:
    if (m_factors.size() == 1) {
        return m_coefficient *
               m_factors[0]->extract_block(row_range, col_range);
    }

    // Iterator over the factors from the back
    auto itfac = m_factors.rbegin();

    // Do the last factor. Note that we will only need those columns
    // desired by our function.
    const auto rows_of_last = range((*itfac)->n_rows());
    stored_matrix_type cache =
          m_coefficient * (*itfac)->extract_block(rows_of_last, col_range);

    // Last factor done:
    ++itfac;

    // Perform the multiplications on the remaining factors one-by-one
    // except the first factor.
    contract_in_place(itfac, m_factors.rend() - 1, cache);

    // Extract the first factor, but only those rows required:
    const auto cols_of_first = range(m_factors.front()->n_cols());
    stored_matrix_type first =
          m_factors.front()->extract_block(row_range, cols_of_first);

    // Return the product of both:
    return first * cache;

    // If extracting the required part of the first factor is much
    // more expansive that doing another lazy*stored multiplication
    // and binning the rows we do not require, than perhaps it is
    // in fact better to let the for loop run over all factors but
    // the last and use another extract_block on the cache to
    // extract the appropriate row elements.
}

template <typename StoredMatrix>
template <typename Iterator>
void LazyMatrixProduct<StoredMatrix>::contract_in_place(
      Iterator begin, Iterator end, stored_matrix_type& m) const {
    for (; begin != end; ++begin) {
        const lazy_matrix_expression_ptr_type& expr_ptr = *begin;
        m = (*expr_ptr) * m;
    }
}

}  // namespace linalgwrap

#endif
