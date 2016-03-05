#pragma once

#include "LazyMatrixExpression.hh"
#include "StoredMatrix_i.hh"
#include <type_traits>
#include "SubscriptionPointer.hh"

namespace linalgwrap {
namespace view {

/** \brief View base interface satisfied by both ViewBase classes*/
template <typename Matrix>
class ViewBase_i {
    typedef Matrix inner_matrix_type;

  public:
    virtual void update(const ParameterMap&) = 0;

  protected:
    /** Access the inner matrix as a reference */
    virtual inner_matrix_type& inner_matrix() = 0;

    /** Access the inner matrix as a const reference */
    virtual const inner_matrix_type& inner_matrix() const = 0;
};

/** \brief View base class for Views of lazy matrices
 *
 * Distinguishing between this class and its brother, the ViewBase for
 * stored matrices is done by the second template argument, which is true_type
 * only for stored matrices.
 * */
template <typename Matrix,
          typename = typename std::is_base_of<
                StoredMatrix_i<typename Matrix::scalar_type>, Matrix>::type>
class ViewBase
      : public LazyMatrixExpression<
              typename std::remove_const<Matrix>::type::stored_matrix_type>,
        public ViewBase_i<Matrix> {
  public:
    typedef Matrix inner_matrix_type;
    typedef LazyMatrixExpression<typename std::remove_const<
          inner_matrix_type>::type::stored_matrix_type> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    static_assert(std::is_base_of<LazyMatrixExpression<stored_matrix_type>,
                                  Matrix>::value,
                  "Matrix is not a child of LazyMatrixExpression of the same "
                  "StoredMatrix type");

    /** A swap function for ViewBase */
    friend void swap(ViewBase& first, ViewBase& second) {
        using std::swap;  // enable ADL

        swap(first.m_inner, second.m_inner);
        swap(first.m_identifier, second.m_identifier);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

    //
    // Constructors, destructors and assignment
    //
    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier An identifier for the type of view.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier)
          : m_inner{inner}, m_identifier(identifier) {}

    ViewBase(const ViewBase&) = default;
    ViewBase(ViewBase&&) = default;

    //
    // Partial implementation of LazyMatrixExpression interface.
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap& map) override {
        CallUpdateIfAvail<Matrix> mat;
        mat(m_inner, map);
        // m_inner.update(map)
    }

    /** \brief Print the expression tree to this outstream */
    void print_tree(std::ostream& o) const override {
        o << m_identifier << " of ";
        inner_matrix().print_tree(o);
    }

  protected:
    /** Access the inner matrix as a const reference */
    const inner_matrix_type& inner_matrix() const override { return m_inner; }

    /** Access the inner matrix as a reference */
    inner_matrix_type& inner_matrix() override { return m_inner; }

  private:
    //! Store a copy of the inner lazy matrix expression:
    inner_matrix_type m_inner;
    std::string m_identifier;
};

/** \brief View base class for Views of stored matrices
 *
 * Distinguishing between this class and its brother, the ViewBase for
 * lazy matrices is done by the second template argument, which is true_type
 * only for stored matrices.
 * */
template <typename Matrix>
class ViewBase<Matrix, std::true_type>
      : public LazyMatrixExpression<typename std::remove_const<Matrix>::type>,
        public ViewBase_i<Matrix> {
  public:
    typedef Matrix inner_matrix_type;
    typedef LazyMatrixExpression<typename std::remove_const<Matrix>::type>
          base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    /** A swap function for ViewBase */
    friend void swap(ViewBase& first, ViewBase& second) {
        using std::swap;  // enable ADL

        swap(first.m_inner_ptr, second.m_inner_ptr);
        swap(static_cast<base_type&>(first), static_cast<base_type&>(second));
    }

    //
    // Constructors, destructors and assignment
    //
    /** Construct a ViewBase class for a LazyMatrix
     *
     * \param inner   Reference to the LazyMatrix
     * \param identifier Identifier for the type of View calling this
     *                   constructor.
     */
    ViewBase(inner_matrix_type& inner, const std::string& identifier)
          : m_inner_ptr{make_subscription(inner, identifier)} {}

    ViewBase(const ViewBase&) = default;
    ViewBase(ViewBase&&) = default;

    //
    // Partial implementation of LazyMatrixExpression interface.
    //
    /** \brief Update the internal data of all objects in this expression
     *         given the ParameterMap
     * */
    void update(const ParameterMap&) override {
        // Do nothing, since the underlying object is a stored
        // matrix.
    }

    /** \brief Print the expression tree to this outstream
     * */
    void print_tree(std::ostream& o) const override {
        o << m_inner_ptr.subscriber_id() << " of " << inner_matrix().name();
    }

  protected:
    /** Access the inner matrix as a reference */
    inner_matrix_type& inner_matrix() override { return *m_inner_ptr; }

    /** Access the inner matrix as a const reference */
    const inner_matrix_type& inner_matrix() const override {
        return *m_inner_ptr;
    }

  private:
    SubscriptionPointer<inner_matrix_type> m_inner_ptr;
};

}  // view
}  // linalgwrap
