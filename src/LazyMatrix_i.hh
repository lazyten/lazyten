#ifndef LINALG_LAZY_MATRIX_I_HPP_
#define LINALG_LAZY_MATRIX_I_HPP_

#include "LazyMatrixExpression.hh"

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

    // Swapping:
    friend void swap(LazyMatrix_i& first, LazyMatrix_i& second) {
        using std::swap;
        swap(first.m_name, second.m_name);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

    //
    // Partial implementation of the interface of a LazyMatrixExpression
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    virtual void update(const ParameterMap&) override {
        assert_dbg(false, ExcNotImplemented());
    }

    void print_tree(std::ostream& o) const override {
        // just print the name of this leaf
        o << name();
    }

    /** Get the name of the matrix */
    std::string name() const { return m_name; }

    /** Set the name of the matrix */
    void name(const std::string& name) { m_name = name; }

  protected:
    //! some name identifying the matrix, or empty
    std::string m_name;
};

}  // namespace linalg
#endif  // LINALG_MATRIX_I_HPP_
