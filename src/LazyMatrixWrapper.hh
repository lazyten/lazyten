#ifndef LINALG_ABSTRACT_MATRIX_HPP_
#define LINALG_ABSTRACT_MATRIX_HPP_

#include "LazyMatrix_i.hh"
#include <memory>

namespace linalgwrap {

/** \brief Class for enabeling lazy evaluation of matrix operations
 *  onto another matrix class, which is not already lazy.
 *
 * \tparam InnerMatrix:   The type of the inner matrix that is made lazy
 * \tparam StoredMatrix:  The type of the stored matrix to use
 */
template <typename StoredMatrix, typename InnerMatrix>
class LazyMatrixWrapper : public LazyMatrixExpression<StoredMatrix> {
    static_assert(
          std::is_same<typename InnerMatrix::scalar_type,
                       typename StoredMatrix::scalar_type>::value,
          "InnerMatrix and StoredMatrix need to have the same scalar type");

    static_assert(
          std::is_base_of<StoredMatrix_i<typename InnerMatrix::scalar_type>,
                          InnerMatrix>::value,
          "InnerMatrix must be a child class of StoredMatrix_i");

  public:
    typedef LazyMatrixExpression<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    typedef InnerMatrix inner_matrix_type;

    //
    // Construction and destruction
    //
    /** \brief Constructor from shared pointer
     *
     * Make a copy of the shared pointer and use that as the inner matrix
     * object. All modifications done via this interface of cause are also
     * noticed from the Pointer provided here.
     */
    explicit LazyMatrixWrapper(std::shared_ptr<inner_matrix_type> inner)
          : m_inner(std::move(inner)) {}

    /** \brief Constructor from inner_matrix_type (for implicit conversion) */
    explicit LazyMatrixWrapper(inner_matrix_type&& inner)
          : m_inner(std::make_shared<inner_matrix_type>(inner)) {}

    /** \brief Constructor from inner_matrix_type (for explicit conversion)
     *
     * Warning: This copies the full content of inner.
     */
    explicit LazyMatrixWrapper(const inner_matrix_type& inner)
          : m_inner(std::make_shared<inner_matrix_type>(inner)) {}

    //
    // Access to inner matrix
    //
    /** Const access to the inner matrix object */
    const inner_matrix_type& inner_matrix() const { return *m_inner; }

    /** Non-const access to the inner data object */
    inner_matrix_type& inner_matrix() { return *m_inner; }

    //
    // Matrix_i interface
    //
    /**
     * See documentation of Matrix_i function of the same name.
     */
    void fill(size_type start_row, size_type start_col,
              SmallMatrix<scalar_type>& block, bool add = false,
              scalar_type c_this = Constants<scalar_type>::one) const override {
        m_inner->fill(start_row, start_col, block, add, c_this);
    }

    /** \brief Number of rows of the matrix */
    size_type n_rows() const override { return m_inner->n_rows(); }

    /** \brief Number of columns of the matrix  */
    size_type n_cols() const override { return m_inner->n_cols(); }

    //
    // LazyMatrixExpression interface
    //
    /** \brief call the update routine of all lazy matrices with the
     *         specified arguments
     */
    template <typename... Args>
    void update(Args...) {
        // Nothing to do here, since static object.
    }

    /** \brief Multiplication with a stored matrix */
    stored_matrix_type operator*(const stored_matrix_type& m) const override {
        // TODO assume that stored_matrix_type and inner_matrix_type
        //      can be multiplied together.
        return (*m_inner) * m;
    }

    /** \brief Print the expression tree to this outstream
     * */
    void print_tree(std::ostream& o) const override {
        // We are the leaf, just print the name of inner:
        o << m_inner->name();
    }

    /** \brief Clone the expression */
    lazy_matrix_expression_ptr_type clone() const override {
        // return a copy enwrapped in the pointer type
        return lazy_matrix_expression_ptr_type(new LazyMatrixWrapper(*this));
    }

  private:
    std::shared_ptr<inner_matrix_type> m_inner;
};

/** \brief Convenience function to construct a LazyMatrixWrapper */
template <typename StoredMatrix, typename InnerMatrix, typename... Args>
LazyMatrixWrapper<StoredMatrix, InnerMatrix> make_lazy_matrix(Args... args) {
    return LazyMatrixWrapper<StoredMatrix, InnerMatrix>(
          std::move(std::make_shared<InnerMatrix>(args...)));
}

}  // namespace linalg

#endif
