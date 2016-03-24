#ifndef LINALG_MATRIX_SUM_HPP_
#define LINALG_MATRIX_SUM_HPP_

#include "Constants.hh"
#include "LazyMatrixExpression.hh"
#include "LazyMatrixProduct.hh"
#include "SubscriptionPointer.hh"
#include <vector>
#include <algorithm>
#include <iterator>

namespace linalgwrap {

// forward declaration
template <typename StoredMatrix>
class LazyMatrixProduct;

template <typename StoredMatrix>
class LazyMatrixExpression;

// TODO: If sparsity should be incorporated, build up the sparsity pattern
//       as more and more terms are added to the sum object.

/** Macro to aid defining both in-place matrix-matrix addition and in-place
 * matrix-matrix subtraction */
#define LazyMatrixSum_inplace_addsub_op(othertype)            \
    LazyMatrixSum& operator+=(othertype other) {              \
        this->push_term(other);                               \
        return *this;                                         \
    }                                                         \
    LazyMatrixSum& operator-=(othertype other) {              \
        this->push_term(other, -Constants<scalar_type>::one); \
        return *this;                                         \
    }

/** Class to represent the sum of different MatrixProducts
 * It may include stored terms as well.
 */
template <typename StoredMatrix>
class LazyMatrixSum : public LazyMatrixExpression<StoredMatrix> {
  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    /** The type of the lazy terms contained in this sum */
    typedef LazyMatrixProduct<stored_matrix_type> lazy_term_type;

  private:
    /** Type used for terms which are stored matrices
     *
     * Stores a reference to a stored matrix (by the means of a
     * subscription pointer) and a coefficient.
     *
     * @note: This implies that changing the stored matrix after
     *        enwrapping it in this term structure changes the
     *        value of this term as well!
     */
    class StoredTerm {
      public:
        StoredTerm(const stored_matrix_type& matrix, scalar_type coefficient)
              : m_coefficient{coefficient},
                m_matrix_ptr{make_subscription(matrix, "LazyMatrixSum")} {}

        //! get the coefficient
        scalar_type coefficient() const { return m_coefficient; }

        //! get the matrix as a reference
        const stored_matrix_type& matrix() const { return *m_matrix_ptr; }

        //! adjust the coefficient
        StoredTerm& operator*=(scalar_type c) {
            m_coefficient *= c;
            return *this;
        }

      private:
        scalar_type m_coefficient;
        SubscriptionPointer<const stored_matrix_type> m_matrix_ptr;
    };

    /** Have an alias for the StoredTerm class */
    typedef StoredTerm stored_term_type;

  public:
    /** A swap function for LazyMatrixSums */
    friend void swap(LazyMatrixSum& first, LazyMatrixSum& second) {
        using std::swap;  // enable ADL

        swap(first.m_n_rows, second.m_n_rows);
        swap(first.m_n_cols, second.m_n_cols);
        swap(first.m_lazy_terms, second.m_lazy_terms);
        swap(first.m_stored_terms, second.m_stored_terms);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

  public:
    /** \name Constructors for lazy sums
     */
    ///@{
    /** \brief Create a matrix sum object
     *
     * @param term   The first matrix expression term
     */
    explicit LazyMatrixSum(lazy_term_type term)
          : m_n_rows(term.n_rows()), m_n_cols(term.n_cols()) {
        m_lazy_terms.push_back(std::move(term));
    }

    /** \brief Create a matrix sum object
     *
     * @param expr    The first matrix expression term
     * @param factor  The factor to multiply the matrix expression with.
     */
    explicit LazyMatrixSum(const LazyMatrixExpression<StoredMatrix>& expr,
                           scalar_type factor = Constants<scalar_type>::one)
          : LazyMatrixSum(lazy_term_type(expr, factor)) {}

    /** \brief Create a matrix sum object
     *
     * In this case a subscription to the term is made, meaning that
     * the term object should not be destroyed before the LazyMatrixSum
     * object. Otherwise an exception will throw in DEBUG mode.
     *
     * @param mat   The stored matrix to use as first term.
     * @param factor The factor to multiply the matrix with
     */
    explicit LazyMatrixSum(const stored_matrix_type& mat,
                           scalar_type factor = Constants<scalar_type>::one)
          : m_n_rows{mat.n_rows()}, m_n_cols{mat.n_cols()} {

        // Construct a term and subscribe to the reference
        stored_term_type term(mat, factor);

        // Push term onto the stored stack
        m_stored_terms.push_back(std::move(term));
    }

    //@{
    /** Defaults for the big five */
    ~LazyMatrixSum() = default;
    LazyMatrixSum(const LazyMatrixSum&) = default;
    LazyMatrixSum(LazyMatrixSum&&) = default;
    LazyMatrixSum& operator=(const LazyMatrixSum&) = default;
    LazyMatrixSum& operator=(LazyMatrixSum&&) = default;
    //@}
    ///@}

    //
    // Push back further terms
    //
    /** Push back a further lazy matrix term */
    void push_term(lazy_term_type term) {
        // Check sizes
        assert_size(n_cols(), term.n_cols());
        assert_size(n_rows(), term.n_rows());

        // Push back
        m_lazy_terms.push_back(std::move(term));
    }

    /** Push back a further lazy matrix expression */
    void push_term(const LazyMatrixExpression<StoredMatrix>& term,
                   scalar_type factor = Constants<scalar_type>::one) {
        // Reduce to version above
        push_term(lazy_term_type(term, factor));
    }

    /** \brief Push back a further stored matrix
     *
     * In this case a subscription to the term is made, meaning that
     * the term object should not be destroyed before the LazyMatrixSum
     * object. Otherwise an exception will throw in DEBUG mode.
     * */
    void push_term(const stored_matrix_type& mat,
                   scalar_type factor = Constants<scalar_type>::one) {
        // Check sizes
        assert_size(n_cols(), mat.n_cols());
        assert_size(n_rows(), mat.n_rows());

        // Construct a term and subscribe to the reference
        stored_term_type term(mat, factor);

        m_stored_terms.push_back(std::move(term));
    }

    void push_term(LazyMatrixSum sum) {
        // Check sizes
        assert_size(n_cols(), sum.n_cols());
        assert_size(n_rows(), sum.n_rows());

        // Move all lazy terms of sum to the end of this object
        std::move(std::begin(sum.m_lazy_terms), std::end(sum.m_lazy_terms),
                  back_inserter(m_lazy_terms));

        // Move all stored terms of sum to the end of this object
        std::move(std::begin(sum.m_stored_terms), std::end(sum.m_stored_terms),
                  back_inserter(m_stored_terms));
    }

    //
    // Other modifications of the sum
    //
    /** \brief scale the whole sum by the scalar \p c */
    void scale(const scalar_type c) {
        assert_finite(c);

        for (auto& term : m_stored_terms) {
            term *= c;
        }

        for (auto& term : m_lazy_terms) {
            term *= c;
        }
    }

    //
    // Matrix_i interface
    //

    /** \brief Number of rows of the matrix      */
    size_type n_rows() const override { return m_n_rows; }

    /** \brief Number of columns of the matrix     */
    size_type n_cols() const override { return m_n_cols; }

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
    stored_matrix_type extract_block(Range<size_type> row_range,
                                     Range<size_type> col_range) const override;

    /** \brief Add a block of values of the matrix to the stored matrix
     *         provided by reference.
     *
     * For more details of the interface see the function of the same name
     * in ``LazyMatrixExpression``.
     *
     *  \param in   Matrix to add the values to. It is assumed that it
     *              already has the correct sparsity structure to take
     *              all the values. Its size defines the size of the
     *              block
     *  \param start_row  The row index of the first element to extract
     *  \param start_col  The column index of the first element to extract
     *  \param c_this     The coefficient to multiply this matrix with
     *                    before extracting.
     */
    void add_block_to(
          stored_matrix_type& in, size_type start_row, size_type start_col,
          scalar_type c_this = Constants<scalar_type>::one) const override;

    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap& map) override {
        // Pass the call onto all factors:
        for (auto& expression : m_lazy_terms) {
            expression.update(map);
        }
    }

    /** \brief Multiplication with a stored matrix */
    virtual stored_matrix_type operator*(
          const stored_matrix_type& m) const override {
        assert_size(n_cols(), m.n_rows());

        // TODO
        // The most common cases with zero or one
        // stored matrix term could perhaps be optimised.

        // Allocate storage for the return and fill with zero
        stored_matrix_type res(n_rows(), m.n_cols(), true);

        for (const auto& stored_term : m_stored_terms) {
            const scalar_type coeff = stored_term.coefficient();
            const stored_matrix_type& mat = stored_term.matrix();

            // do the multiplication:
            const stored_matrix_type matm = mat * m;

            // res += coeff*matm
            matm.add_block_to(res, 0, 0, coeff);
        }

        for (const auto& term : m_lazy_terms) {
            // evaluate term times m and add to result:
            res += term * m;
        }
        return res;
    }

    /** \brief Print the expression tree to this outstream
     * */
    virtual void print_tree(std::ostream& o) const override {
        // TODO to be implemented
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixSum(*this));
    }

    //
    // In-place scaling operators
    //
    /** \brief scale a matrix using a scalar
     */
    LazyMatrixSum& operator*=(scalar_type c) {
        assert_finite(c);
        this->scale(c);
        return *this;
    }

    /** \brief scale a matrix using a scalar
     */
    LazyMatrixSum& operator/=(scalar_type c) {
        assert_finite(c);
        assert_dbg(c != 0, ExcDevideByZero());
        this->scale(Constants<scalar_type>::one / c);
        return *this;
    }

    //
    // Operators with += or -=
    //
    // TODO rewrite macros such that clang and YCM are both happy
    //      and do not show random errors
    LazyMatrixSum_inplace_addsub_op(lazy_term_type);
    LazyMatrixSum_inplace_addsub_op(const LazyMatrixExpression<StoredMatrix>&);
    LazyMatrixSum_inplace_addsub_op(const stored_matrix_type&);
    LazyMatrixSum_inplace_addsub_op(LazyMatrixSum);

  private:
    //! The cached number of rows
    size_type m_n_rows;

    //! The cached number of columns
    size_type m_n_cols;

    //! Collection of all the terms which are subject to delayed evaluation
    std::vector<lazy_term_type> m_lazy_terms;

    /** Collection of stored matrix terms to which we reference */
    std::vector<stored_term_type> m_stored_terms;
};

//
// Operators which perform scaling
//
/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator*(LazyMatrixSum<StoredMatrix> lhs,
                                      typename StoredMatrix::scalar_type rhs) {
    lhs *= rhs;
    return lhs;
}

/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator*(typename StoredMatrix::scalar_type lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
    return rhs * lhs;
}

/** \brief Scale a lazy matrix sum */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator/(LazyMatrixSum<StoredMatrix> lhs,
                                      typename StoredMatrix::scalar_type rhs) {
    lhs /= rhs;
    return lhs;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> mat) {
    typedef typename StoredMatrix::scalar_type scalar_type;
    return -Constants<scalar_type>::one * mat;
}

//
// Macro definition for adding and subtracting things to/from LazyMatrixSums
//
#define LazyMatrixSum_outofplace_addsub_op(othertype)                        \
    template <typename StoredMatrix>                                         \
    LazyMatrixSum<StoredMatrix> operator+(LazyMatrixSum<StoredMatrix> sum,   \
                                          othertype other) {                 \
        sum += other;                                                        \
        return sum;                                                          \
    }                                                                        \
    template <typename StoredMatrix>                                         \
    LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> sum,   \
                                          othertype other) {                 \
        sum -= other;                                                        \
        return sum;                                                          \
    }                                                                        \
    template <typename StoredMatrix>                                         \
    LazyMatrixSum<StoredMatrix> operator+(othertype other,                   \
                                          LazyMatrixSum<StoredMatrix> sum) { \
        return sum + other;                                                  \
    }                                                                        \
    template <typename StoredMatrix>                                         \
    LazyMatrixSum<StoredMatrix> operator-(othertype other,                   \
                                          LazyMatrixSum<StoredMatrix> sum) { \
        return (-sum) + other;                                               \
    }

//
// Macro definition for adding and subtracting things to/from LazyMatrixProducts
//
#define LazyMatrixProduct_outofplace_addsub_op(othertype)          \
    template <typename StoredMatrix>                               \
    LazyMatrixSum<StoredMatrix> operator+(                         \
          LazyMatrixProduct<StoredMatrix> prod, othertype other) { \
        LazyMatrixSum<StoredMatrix> sum(std::move(prod));          \
        return sum + other;                                        \
    }                                                              \
    template <typename StoredMatrix>                               \
    LazyMatrixSum<StoredMatrix> operator-(                         \
          LazyMatrixProduct<StoredMatrix> prod, othertype other) { \
        LazyMatrixSum<StoredMatrix> sum(std::move(prod));          \
        return sum - other;                                        \
    }                                                              \
    template <typename StoredMatrix>                               \
    LazyMatrixSum<StoredMatrix> operator+(                         \
          othertype other, LazyMatrixProduct<StoredMatrix> prod) { \
        return prod + other;                                       \
    }                                                              \
    template <typename StoredMatrix>                               \
    LazyMatrixSum<StoredMatrix> operator-(                         \
          othertype other, LazyMatrixProduct<StoredMatrix> prod) { \
        return (-prod) + other;                                    \
    }

//
// Further addition operators
//
LazyMatrixSum_outofplace_addsub_op(LazyMatrixProduct<StoredMatrix>);
LazyMatrixSum_outofplace_addsub_op(const StoredMatrix&);
LazyMatrixSum_outofplace_addsub_op(const LazyMatrixExpression<StoredMatrix>&);

LazyMatrixProduct_outofplace_addsub_op(const StoredMatrix&);
LazyMatrixProduct_outofplace_addsub_op(
      const LazyMatrixExpression<StoredMatrix>&);

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(LazyMatrixSum<StoredMatrix> lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
    lhs += rhs;
    return lhs;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixSum<StoredMatrix> lhs,
                                      LazyMatrixSum<StoredMatrix> rhs) {
    lhs -= rhs;
    return lhs;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(LazyMatrixProduct<StoredMatrix> lhs,
                                      LazyMatrixProduct<StoredMatrix> rhs) {
    LazyMatrixSum<StoredMatrix> sum(std::move(lhs));
    sum += rhs;
    return sum;
}

template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(LazyMatrixProduct<StoredMatrix> lhs,
                                      LazyMatrixProduct<StoredMatrix> rhs) {
    LazyMatrixSum<StoredMatrix> sum(std::move(lhs));
    sum -= rhs;
    return sum;
}

//
// ------------------------------------------------------------
//
//
// LazyMatrixSum
//
template <typename StoredMatrix>
typename LazyMatrixSum<StoredMatrix>::scalar_type LazyMatrixSum<StoredMatrix>::
operator()(size_type row, size_type col) const {
    assert_range(0, row, n_rows());
    assert_range(0, col, n_cols());
    auto block = extract_block({row, row + 1}, {col, col + 1});
    return block(0, 0);
}

template <typename StoredMatrix>
typename LazyMatrixSum<StoredMatrix>::stored_matrix_type
LazyMatrixSum<StoredMatrix>::extract_block(Range<size_type> row_range,
                                           Range<size_type> col_range) const {
    // At least one range is empty -> no work to be done:
    if (row_range.is_empty() || col_range.is_empty()) {
        return stored_matrix_type{row_range.length(), col_range.length()};
    }

    // Assertive checks:
    assert_greater_equal(row_range.last(), this->n_rows());
    assert_greater_equal(col_range.last(), this->n_cols());

    // Allocate storage and set elements to zero
    stored_matrix_type res(row_range.length(), col_range.length(), true);

    // Add all terms to res
    add_block_to(res, row_range.first(), col_range.first());

    // Return it
    return res;
}

template <typename StoredMatrix>
void LazyMatrixSum<StoredMatrix>::add_block_to(stored_matrix_type& in,
                                               size_type start_row,
                                               size_type start_col,
                                               scalar_type c_this) const {
    // check that we do not overshoot the row index
    assert_greater_equal(start_row + in.n_rows(), this->n_rows());

    // check that we do not overshoot the column index
    assert_greater_equal(start_col + in.n_cols(), this->n_cols());

    for (const auto& stored_term : m_stored_terms) {
        const scalar_type coeff = stored_term.coefficient();
        const stored_matrix_type& mat = stored_term.matrix();

        // Add the elements from the scaled matrix to the present
        // block
        mat.add_block_to(in, start_row, start_col, c_this * coeff);
    }

    for (const auto& term : m_lazy_terms) {
        // add results to the present block:
        term.add_block_to(in, start_row, start_col, c_this);
    }
}
}  // namespace linalg
#endif
