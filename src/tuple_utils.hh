#pragma once
#include "Constants.hh"
#include <functional>
#include <tuple>
#include <utility>

// TODO try to make as much use of the new c++14 index_sequence object here
//      and get rid of the head/tail stuff as much as possible, as it makes
//      it harder for the compiler to optimise these operations.

namespace linalgwrap {
namespace detail {
//! Implementation to get a tuple of references to the tail of a tuple.
template <std::size_t... Indices, typename E, typename... Es>
std::tuple<Es&...> tuple_tail_impl(std::index_sequence<Indices...>,
                                   std::tuple<E&, Es&...> t);

//! Implementation of the ref function for tuple references
//  We want to return a tuple of references from a referenced tuple.
template <std::size_t... Indices, typename... Es>
std::tuple<Es&...> tuple_ref_impl(std::index_sequence<Indices...>,
                                  std::tuple<Es...>& t);
}  // namespace detail

/** \name Basic tuple operations */
///@{
/** Get the head of a std::tuple of references
 *
 * \note You can use tuple_ref to make a tuple of references out of a
 *       referenced tuple
 * */
template <typename E, typename... Es>
E& tuple_head(std::tuple<E&, Es&...> t);

/** Get the tail of a std::tuple of references
 *
 * \note You can use tuple_ref to make a tuple of references out of a
 *       referenced tuple
 * */
template <typename E, typename... Es>
std::tuple<Es&...> tuple_tail(std::tuple<E&, Es&...> t);

/** Construct a tuple taking ownership of head and tail. */
template <typename E, typename... Es>
std::tuple<E, Es...> tuple_cons(E&& e, std::tuple<Es...>&& t);

/** Construct a tuple of references to each element of the referenced tuple */
template <typename... Es>
std::tuple<Es&...> tuple_ref(std::tuple<Es...>& t);
///@}

/** \name Tuple traversal and mapping */
///@{
//@{
/** \brief Apply the operation op to the first element
 *         in a tuple of references, which matches a predicate.
 *
 *         Does not apply the operation if no object
 *         matches the predicate
 *
 * \note You can use tuple_ref to make a tuple of references out of a
 *       referenced tuple
 */
template <typename Pred, typename Op, typename... Ts>
void tuple_for_first(Pred pred, Op op, std::tuple<Ts&...> t);

// End recursion
template <typename Pred, typename Op>
void tuple_for_first(Pred, Op, std::tuple<>);
//@}

//@{
/** \brief Apply a unary operator to all elements of a tuple of references
 *         in turn.
 *
 * \note You can use tuple_ref to make a tuple of references out of a
 *       referenced tuple
 **/
template <typename UnOp, typename... Ts>
void tuple_for_each(UnOp op, std::tuple<Ts&...> t);

// Special case to end recursion
template <typename UnOp>
void tuple_for_each(UnOp, std::tuple<>);
//@}

//@{
/** \brief Map a unary operator to all elements of a tuples of references
 *         in turn.
 *
 * Return the tuple of all results.
 *
 *  Of cause the apply function of the operator should be generic in
 *  all types occurring.
 *
 * \note The operator has to return a non-void object, otherwise this
 * implementation fails. Use tuple_for_each for unary operators
 * returning void.
 *
 * \note You can use tuple_ref to make a tuple of references out of a
 *       referenced tuple
 *  */
template <typename UnOp, typename... Ts>
auto tuple_map(UnOp op, std::tuple<Ts&...> t);

// Special case to end recursion.
template <typename UnOp>
std::tuple<> tuple_map(UnOp, std::tuple<>);
//@}

//@{
/** Map a binary operator to all elements of two tuples in turn.
 *
 * Return the tuple of all results.
 *
 *  Of cause the apply function of the operator should be generic in
 *  all pairs of types occurring.
 *
 * \note The operator has to return a non-void object, otherwise this
 * implementation fails. Use for_each for operators
 * returning void.
 *  */
// Special case to end recursion.
template <typename BinOp>
std::tuple<> tuple_map(BinOp, std::tuple<>, std::tuple<>);

template <typename BinOp, typename... Tlhs, typename... Trhs>
auto tuple_map(BinOp op, std::tuple<Tlhs&...> lhs, std::tuple<Trhs&...> rhs);
//@}
///@}

//
// ---------------------------------------------------
//

namespace detail {
template <std::size_t... Indices, typename E, typename... Es>
inline std::tuple<Es&...> tuple_tail_impl(std::index_sequence<Indices...>,
                                          std::tuple<E&, Es&...> t) {
    // Return a tuple of the references we extract with get.
    return std::tie(std::get<Indices + 1u>(t)...);
}

template <std::size_t... Indices, typename... Es>
inline std::tuple<Es&...> tuple_ref_impl(std::index_sequence<Indices...>,
                                         std::tuple<Es...>& t) {
    // Return a tuple of the references we extract with get.
    return std::tie(std::get<Indices>(t)...);
}
}  // namespace detail

template <typename E, typename... Es>
inline E& tuple_head(std::tuple<E&, Es&...> t) {
    return std::get<0>(t);
}

/** Get the tail of a std::tuple reference*/
template <typename E, typename... Es>
inline std::tuple<Es&...> tuple_tail(std::tuple<E&, Es&...> t) {
    // First construct an index sequence from 0 to the number of types in Es
    auto idcs = std::make_index_sequence<sizeof...(Es)>();

    // Then call the implementation function
    return detail::tuple_tail_impl(idcs, t);
}

template <typename E, typename... Es>
inline std::tuple<E, Es...> tuple_cons(E&& e, std::tuple<Es...>&& t) {
    return std::tuple_cat(std::tuple<E>{e}, t);
}

template <typename... Es>
std::tuple<Es&...> tuple_ref(std::tuple<Es...>& t) {
    auto idcs = std::make_index_sequence<sizeof...(Es)>();
    return detail::tuple_ref_impl(idcs, t);
}

template <typename Pred, typename Op, typename... Ts>
inline void tuple_for_first(Pred pred, Op op, std::tuple<Ts&...> t) {
    if (pred(tuple_head(t))) {
        op(tuple_head(t));
        return;
    } else {
        tuple_for_first(std::forward<Pred>(pred), std::forward<Op>(op), t);
    }
}

template <typename Pred, typename Op>
inline void tuple_for_first(Pred, Op, std::tuple<>) {}

template <typename UnOp, typename... Ts>
inline void tuple_for_each(UnOp op, std::tuple<Ts&...> t) {
    // Apply to head:
    op(tuple_head(t));

    // Do operation on tail:
    return tuple_for_each(std::forward<UnOp>(op), tuple_tail(t));
}

template <typename UnOp>
inline void tuple_for_each(UnOp, std::tuple<>) {}

template <typename UnOp, typename... Ts>
inline auto tuple_map(UnOp op, std::tuple<Ts&...> t) {
    // Apply op to head:
    auto res_head = op(head(t));

    // Do operation (recursively) on the tail:
    auto res_tail = tuple_map(std::forward<UnOp>(op), tail(t));

    // Cons the results together:
    return tuple_cons(std::move(res_head), std::move(res_tail));
}

// Special case to end recursion.
template <typename UnOp>
inline std::tuple<> tuple_map(UnOp, std::tuple<>) {
    return std::tuple<>{};
}

template <typename BinOp>
inline std::tuple<> tuple_map(BinOp, std::tuple<>, std::tuple<>) {
    return std::tuple<>{};
}

template <typename BinOp, typename... Tlhs, typename... Trhs>
inline auto tuple_map(BinOp op, std::tuple<Tlhs&...> lhs,
                      std::tuple<Trhs&...> rhs) {
    static_assert(sizeof...(Tlhs) == sizeof...(Trhs),
                  "Both tuples to map for the case of binary operations "
                  "need to be of the same size.");

    // Extract head as references
    auto& lhs_head = tuple_head(lhs);
    auto& rhs_head = tuple_head(rhs);

    // Do the operation on the head:
    auto res_head = op(lhs_head, rhs_head);

    // Extract tail (as tuple of references)
    auto lhs_tail = tuple_tail(lhs);
    auto rhs_tail = tuple_tail(rhs);

    // Do operation (recursively) on the tail:
    auto res_tail = tuple_map(std::forward<BinOp>(op), lhs_tail, rhs_tail);

    // Cons the results together:
    return tuple_cons(std::move(res_head), std::move(res_tail));
}

}  // namespace linalgwrap
