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
 *
 * TODO implement optimisation if two stored matrices are to be multiplied
 *      directly after another (should not happen automatically, so not of
 *      major concern).
 *
 * TODO One could do this by keeping track via a speciffic pointer whether
 *      the last object pushed onto the factor stack was in fact a stored
 *      matrix enwraped into a lazy matrix and if yes do the multiplication
 *      straight away.
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

    //
    // Constructors, destructors, assignment
    //
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
    void scale(const scalar_type c) { m_coefficient *= c; }

    //
    // Matrix_i interface
    //
    /**
     * See documentation of Matrix_i function of the same name.
     */
    // TODO comment: this is expensive for products of more than 2 lazy matrices
    virtual void fill(size_type start_row, size_type start_col,
                      SmallMatrix<scalar_type>& block, bool add = false,
                      scalar_type c_this = Constants<scalar_type>::one) const {
        // check that we do not overshoot the row index
        assert_upper_bound(start_row + block.n_rows() - 1, n_rows());

        // check that we do not overshoot the column index
        assert_upper_bound(start_col + block.n_cols() - 1, n_cols());

        // if we have only one factor, just do it there and be done
        if (m_factors.size() == 1) {
            // Pass on to the factor, just incorporating the c_this
            // and the m_coefficient.
            m_factors[0]->fill(start_row, start_col, block, add,
                               m_coefficient * c_this);
            return;
        }

        auto it = std::begin(m_factors);

        // From the first factor get the required rows from all columns
        SmallMatrix<scalar_type> cache(block.n_rows(), (*it)->n_cols(), false);
        (*it)->fill(start_row, 0, cache);

        // first factor done:
        ++it;

        for (; it < (std::end(m_factors) - 1); ++it) {
            // Get all rows and contract them away
            SmallMatrix<scalar_type> currentmat((*it)->n_rows(),
                                                (*it)->n_cols(), false);
            (*it)->fill(0, 0, currentmat);
            cache = cache * currentmat;
        }

        // From the last factor get the required cols from all rows
        SmallMatrix<scalar_type> lastcache((*it)->n_rows(), block.n_cols(),
                                           false);
        (*it)->fill(0, start_col, lastcache);

        // finally calculate the result.
        if (add) {
            block += m_coefficient * c_this * cache * lastcache;
        } else {
            block = m_coefficient * c_this * cache * lastcache;
        }
    }

    /** \brief Number of rows of the matrix */
    virtual size_type n_rows() const {
        // The number of rows of the product
        // equals the number of rows of the first element in the product.
        return m_factors.front()->n_rows();
    }

    /** \brief Number of columns of the matrix  */
    virtual size_type n_cols() const {
        // The number of columns of the product
        // equals the nuber of columns of the last element in the product.
        return m_factors.back()->n_cols();
    }

    //
    // LazyMatrixExpression interface
    //
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
    virtual stored_matrix_type operator*(const stored_matrix_type& m) const {
        assert_size(n_cols(), m.n_rows());

        // Deal with last factor:
        auto it = m_factors.rbegin();
        stored_matrix_type res = m_coefficient * (*(*it) * m);
        ++it;

        // Deal with all others:
        for (; it != m_factors.rend(); ++it) {
            // Pointer to the factor object
            const lazy_matrix_expression_ptr_type& expr_ptr = *it;
            res = (*expr_ptr) * res;
        }
        return res;
    }

    /** \brief Print the expression tree to this outstream
     * */
    virtual void print_tree(std::ostream& o) const {
        // TODO to be implemented
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixProduct(*this));
    }

    //
    // In-place scalar operators:
    //
    /** \brief scale a matrix using a scalar
     */
    LazyMatrixProduct& operator*=(scalar_type c) {
        this->scale(c);
        return *this;
    }

    LazyMatrixProduct& operator/=(scalar_type c) {
        this->scale(Constants<scalar_type>::one / c);
        return *this;
    }

  private:
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

}  // namespace linalgwrap

#endif
