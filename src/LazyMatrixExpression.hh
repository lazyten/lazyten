#ifndef LINALG_LAZY_MATRIX_EXPRESSION_HPP_
#define LINALG_LAZY_MATRIX_EXPRESSION_HPP_

#include "StoredMatrix_i.hh"
#include "Matrix_i.hh"
#include "LazyMatrixSum.hh"
#include "LazyMatrixProduct.hh"
#include <ostream>
#include <memory>

namespace linalgwrap {

// Forward declarations
template <typename StoredMatrix>
class LazyMatrixSum;

template <typename StoredMatrix>
class LazyMatrixProduct;

/** \brief Generic LazyMatrixExpression class
 *
 * Abstract base class for all lazy matrix expressions.
 * Implements a slightly modified paradigm of Expression templates.
 *
 * \tparam StoredMatrix The type to use for stored matrices.
 * */
template <typename StoredMatrix>
class LazyMatrixExpression
      : public Matrix_i<typename StoredMatrix::scalar_type> {
    static_assert(
          std::is_base_of<StoredMatrix_i<typename StoredMatrix::scalar_type>,
                          StoredMatrix>::value,
          "StoredMatrix is not a child of StoredMatrix_i of the same scalar "
          "type");

  public:
    typedef StoredMatrix stored_matrix_type;
    typedef Matrix_i<typename StoredMatrix::scalar_type> base_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::size_type size_type;

    /** The pointer type to use for pointers to LazyMatrixExpressions */
    typedef std::shared_ptr<LazyMatrixExpression<StoredMatrix>>
          lazy_matrix_expression_ptr_type;

    // Swapping:
    friend void swap(LazyMatrixExpression&, LazyMatrixExpression&) {
        // nothing
    }

    //
    // Extra interface for expressions:
    //
    /** \brief Print the expression tree to this outstream */
    virtual void print_tree(std::ostream& o) const = 0;

    /** \brief Update the internal data of all objects in this expression
     * given a set of args
     *
     * TODO Cannot be a template since it needs to be virtual
     * */
    template <typename... Args>
    void update(Args...) {
        assert_dbg(false, ExcNotImplemented());
    }

    /** \brief Multiplication with a stored matrix */
    virtual stored_matrix_type operator*(
          const stored_matrix_type& in) const = 0;

    /** \brief Clone the expression
     *
     * Return a clone of the current object as a pointer of type
     * lazy_matrix_expression_ptr_type
     */
    virtual lazy_matrix_expression_ptr_type clone() const = 0;
};

//
// Multiplication
//
/** Multiply two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct product with one factor:
    LazyMatrixProduct<StoredMatrix> prod{lhs};

    // add the other factor and return the result.
    prod.push_factor(rhs);
    return prod;
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      const LazyMatrixExpression<StoredMatrix>& m,
      typename StoredMatrix::scalar_type s) {

    // Construct product with one factor:
    LazyMatrixProduct<StoredMatrix> prod{m};

    // and scale it:
    prod *= s;
    return prod;
}

/** Scale a lazy matrix expression */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator*(
      typename StoredMatrix::scalar_type s,
      const LazyMatrixExpression<StoredMatrix>& m) {
    return m * s;
}

//
// Division by scalar
//
/** \brief Devide a lazy matrix by a scalar  */
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator/(
      const LazyMatrixExpression<StoredMatrix>& m,
      typename StoredMatrix::scalar_type s) {

    typedef typename StoredMatrix::scalar_type scalar_type;
    scalar_type inverse = Constants<scalar_type>::one / s;
    return m * inverse;
}

//
// Addition
//
/** Add two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // add the other factor and return the result.
    sum.push_term(rhs);
    return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const LazyMatrixExpression<StoredMatrix>& lhs, const StoredMatrix& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // add the other factor and return the result.
    sum.push_term(rhs);
    return sum;
}

/** Add a lazy matrix and a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator+(
      const StoredMatrix& lhs, const LazyMatrixExpression<StoredMatrix>& rhs) {
    return rhs + lhs;
}

//
// Subtraction
//
/** Subtract two lazy matrix expressions */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& lhs,
      const LazyMatrixExpression<StoredMatrix>& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // subtract the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

/** Subtract a stored matrix from a lazy matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& lhs, const StoredMatrix& rhs) {

    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // add the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

/** Subtract a lazy matrix from a stored matrix */
template <typename StoredMatrix>
LazyMatrixSum<StoredMatrix> operator-(
      const StoredMatrix& lhs, const LazyMatrixExpression<StoredMatrix>& rhs) {
    // Construct sum with one term:
    LazyMatrixSum<StoredMatrix> sum{lhs};

    // typedef the scalar type
    typedef typename StoredMatrix::scalar_type scalar_type;

    // add the other factor and return the result.
    sum.push_term(rhs, -Constants<scalar_type>::one);
    return sum;
}

//
// Unary operator-
//
template <typename StoredMatrix>
LazyMatrixProduct<StoredMatrix> operator-(
      const LazyMatrixExpression<StoredMatrix>& mat) {
    typedef typename StoredMatrix::scalar_type scalar_type;
    return -Constants<scalar_type>::one * mat;
}

}  // namespace linalg
#endif  // LINALG_LAZY_MATRIX_EXPRESSION_HPP_
