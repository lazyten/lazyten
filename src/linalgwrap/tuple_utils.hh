#pragma once
#include "linalgwrap/Constants.hh"
#include <functional>
#include <tuple>
#include <utility>

// TODO try to make as much use of the new c++14 index_sequence object here
//      and get rid of the head/tail stuff as much as possible, as it makes
//      it harder for the compiler to optimise these operations.

namespace linalgwrap {
namespace detail {
//! Implementation to get the tuple tail
template <std::size_t... Indices, typename E, typename... Es>
constexpr std::tuple<Es...> tuple_tail_impl(std::index_sequence<Indices...>,
                                            std::tuple<E, Es...>&& t);

//! Implementation of the ref function for tuple references
//  We want to return a tuple of references from a referenced tuple.
template <std::size_t... Indices, typename... Es>
constexpr std::tuple<Es&...> tuple_ref_impl(std::index_sequence<Indices...>,
                                            std::tuple<Es...>& t);
}  // namespace detail

/** \name Basic tuple operations */
///@{
/** Get the head of a std::tuple. */
template <typename E, typename... Es>
constexpr E& tuple_head(std::tuple<E, Es...>& t);

/** Get the head of a const std::tuple. */
template <typename E, typename... Es>
constexpr const E& tuple_head(const std::tuple<E, Es...>& t);

/** Get the tail of a std::tuple, taking ownership of the original object. */
template <typename E, typename... Es>
constexpr std::tuple<Es...> tuple_tail(std::tuple<E, Es...>&& t);

/** Construct a tuple taking ownership of head and tail. */
template <typename E, typename... Es>
constexpr std::tuple<E, Es...> tuple_cons(E&& e, std::tuple<Es...>&& t);

/** Construct a tuple of references to each element of the referenced tuple */
template <typename... Es>
constexpr std::tuple<Es&...> tuple_ref(std::tuple<Es...>& t);
///@}

/** \name Tuple traversal and mapping */
///@{
//@{
/** \brief Apply the operation op to the first element
 *         in a tuple, which matches a predicate.
 *
 *         Does not apply the operation if no object
 *         matches the predicate
 */
template <typename Pred, typename Op, typename... Ts>
constexpr void tuple_for_first(Pred&& pred, Op&& op, std::tuple<Ts...>&& t);

// End recursion
template <typename Pred, typename Op>
constexpr void tuple_for_first(Pred&&, Op&&, std::tuple<>);
//@}

//@{
/** \brief Apply a unary operator to all elements of a tuple in turn. */
template <typename UnOp, typename... Ts>
constexpr void tuple_for_each(UnOp&& op, std::tuple<Ts...>&& t);

// Special case to end recursion
template <typename UnOp>
constexpr void tuple_for_each(UnOp&&, std::tuple<>);
//@}

//@{
/** \brief Map a unary operator to all elements of a tuple  in turn.
 *
 * Return the tuple of all results.
 *
 *  Of cause the operator() of the UnOp should be generic in
 *  all types occurring.
 *
 * \note The operator has to return a non-void object, otherwise this
 * implementation fails. Use tuple_for_each for unary operators
 * returning void.
 *  */
template <typename UnOp, typename... Ts>
constexpr auto tuple_map(UnOp&& op, std::tuple<Ts...>&& t);

// Special case to end recursion.
template <typename UnOp>
constexpr std::tuple<> tuple_map(UnOp&&, std::tuple<>);
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
template <typename BinOp, typename... Tlhs, typename... Trhs>
constexpr auto tuple_map(BinOp&& op, std::tuple<Tlhs...>&& lhs,
                         std::tuple<Trhs...>&& rhs);

// Special case to end recursion.
template <typename BinOp>
constexpr std::tuple<> tuple_map(BinOp&&, std::tuple<>, std::tuple<>);
//@}
///@}

//
// ---------------------------------------------------
//

namespace detail {
template <std::size_t... Indices, typename E, typename... Es>
constexpr std::tuple<Es...> tuple_tail_impl(std::index_sequence<Indices...>,
                                            std::tuple<E, Es...>&& t) {
    typedef std::tuple<E, Es...> tuple_t;
    typedef std::tuple<Es...> tail_t;

    // Return a tuple with the first element removed
    return tail_t{std::get<Indices + 1u>(std::forward<tuple_t>(t))...};
}

template <std::size_t... Indices, typename... Es>
constexpr std::tuple<Es&...> tuple_ref_impl(std::index_sequence<Indices...>,
                                            std::tuple<Es...>& t) {
    // Return a tuple of the references we extract with get.
    return std::tie(std::get<Indices>(t)...);
}
}  // namespace detail

template <typename E, typename... Es>
constexpr E& tuple_head(std::tuple<E, Es...>& t) {
    return std::get<0>(t);
}

template <typename E, typename... Es>
constexpr const E& tuple_head(const std::tuple<E, Es...>& t) {
    return std::get<0>(t);
}

/** Get the tail of a std::tuple reference*/
template <typename E, typename... Es>
constexpr std::tuple<Es...> tuple_tail(std::tuple<E, Es...>&& t) {
    typedef std::tuple<E, Es...> tuple_t;

    // First construct an index sequence from 0 to the number of types in Es
    auto idcs = std::make_index_sequence<sizeof...(Es)>();

    // Then call the implementation function
    return detail::tuple_tail_impl(idcs, std::forward<tuple_t>(t));
}

template <typename E, typename... Es>
constexpr std::tuple<E, Es...> tuple_cons(E&& e, std::tuple<Es...>&& t) {
    return std::tuple_cat(std::tuple<E>{e}, t);
}

template <typename... Es>
constexpr std::tuple<Es&...> tuple_ref(std::tuple<Es...>& t) {
    auto idcs = std::make_index_sequence<sizeof...(Es)>();
    return detail::tuple_ref_impl(idcs, t);
}

template <typename Pred, typename Op, typename... Ts>
constexpr void tuple_for_first(Pred&& pred, Op&& op, std::tuple<Ts...>&& t) {
    typedef std::tuple<Ts...> tuple_t;

    if (pred(tuple_head(t))) {
        op(tuple_head(t));
        return;
    } else {
        // Extract tail:
        auto tail = tuple_tail(std::forward<tuple_t>(t));

        // Recurse into tail:
        tuple_for_first(std::forward<Pred>(pred), std::forward<Op>(op),
                        std::move(tail));
    }
}

template <typename Pred, typename Op>
constexpr void tuple_for_first(Pred&&, Op&&, std::tuple<>) {}

template <typename UnOp, typename... Ts>
constexpr void tuple_for_each(UnOp&& op, std::tuple<Ts...>&& t) {
    typedef std::tuple<Ts...> tuple_t;

    // Apply to head:
    op(tuple_head(t));

    // Extract tail:
    auto tail = tuple_tail(std::forward<tuple_t>(t));

    // Do operation on tail:
    return tuple_for_each(std::forward<UnOp>(op), std::move(tail));
}

template <typename UnOp>
constexpr void tuple_for_each(UnOp&&, std::tuple<>) {}

template <typename UnOp, typename... Ts>
constexpr auto tuple_map(UnOp&& op, std::tuple<Ts...>&& t) {
    typedef std::tuple<Ts...> tuple_t;

    // Apply op to head:
    auto res_head = op(head(t));

    // Extract tail:
    auto tail = tuple_tail(std::forward<tuple_t>(t));

    // Do operation (recursively) on the tail:
    auto res_tail = tuple_map(std::forward<UnOp>(op), std::move(tail));

    // Cons the results together:
    return tuple_cons(std::move(res_head), std::move(res_tail));
}

// Special case to end recursion.
template <typename UnOp>
constexpr std::tuple<> tuple_map(UnOp&&, std::tuple<>) {
    return std::tuple<>{};
}

template <typename BinOp, typename... Tlhs, typename... Trhs>
constexpr auto tuple_map(BinOp&& op, std::tuple<Tlhs...>&& lhs,
                         std::tuple<Trhs...>&& rhs) {
    static_assert(sizeof...(Tlhs) == sizeof...(Trhs),
                  "Both tuples to map for the case of binary operations "
                  "need to be of the same size.");

    // Extract head as references
    auto& lhs_head = tuple_head(lhs);
    auto& rhs_head = tuple_head(rhs);

    // Do the operation on the head:
    auto res_head = op(lhs_head, rhs_head);

    // Extract tail (as tuple of references)
    auto lhs_tail = tuple_tail(std::forward<std::tuple<Tlhs...>>(lhs));
    auto rhs_tail = tuple_tail(std::forward<std::tuple<Trhs...>>(rhs));

    // Do operation (recursively) on the tail:
    auto res_tail = tuple_map(std::forward<BinOp>(op), std::move(lhs_tail),
                              std::move(rhs_tail));

    // Cons the results together:
    return tuple_cons(std::move(res_head), std::move(res_tail));
}

template <typename BinOp>
constexpr std::tuple<> tuple_map(BinOp&&, std::tuple<>, std::tuple<>) {
    return std::tuple<>{};
}

}  // namespace linalgwrap
