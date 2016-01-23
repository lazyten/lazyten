#ifndef SUBSCRIBABLE_H_
#define SUBSCRIBABLE_H_

#include "Exceptions.hh"
#include "SubscriptionPointer.hh"
#include <string>
#include <list>
#include <typeinfo>
#include <utility>
#include <sstream>

namespace linalgwrap {

// forwand declare SubscriptionPointer
template <typename Subscribable>
class SubscriptionPointer;

/** Class which handles subscribtions from a subscription pointer
 * If upon deletion still subscriptions exist, it throws an exception in debug
 * mode
 */
class Subscribable {
    // Declare SubscriptionPointer as friend.
    template <typename T>
    friend class ::linalgwrap::SubscriptionPointer;

  public:
#ifdef DEBUG
    /** Exception to indicate that Subscribable is still used */
    DefException3(ExcStillUsed, std::string&, size_t, std::string&,
                  << "Object of type " << arg1 << " is still used by " << arg2
                  << " other objects, which are (from new to old): " << arg3);

    /** Exception to indicate that Subscriber is not known. */
    DefException2(ExcUnknownSubscriberId, std::string&, std::string&,
                  << "No subscriber with identifier " << arg1
                  << " is known to have subscribed to the class " << arg2
                  << ".");
#endif

    /** Check if all subscriptions have been removed
     *
     * @note This function is empty unless we are in DEBUG mode
     */
    virtual ~Subscribable() {
#ifdef DEBUG
        if (m_subscribers.size() > 0) {
            // build the string of subscribing objects
            std::stringstream s;
            for (auto& p : m_subscribers) {
                s << " " << p.second;
            }

            // Raise the exception
            assert_dbg(false, ExcStillUsed(m_classname, m_subscribers.size(),
                                           s.str()));
        }
#endif
    }

    /** Return the current number of subscriptions to this object.
     *
     * @note If we are not in DEBUG mode this count is always zero.
     * */
    size_t n_subscriptions() const {
#ifdef DEBUG
        // Note: Getting the size of a list in C++11 is constant time!
        return m_subscribers.size();
#else
        // return constant 0
        return 0;
#endif
    }

  private:
    /** Remove a subscription.
     *
     * @param id Reference to the same string object which was used upon
     *           subscription
     *
     * @note only has an effect if we are in DEBUG mode
     * */
    void unsubscribe(const std::string& id) const {
#ifdef DEBUG
        for (auto it = m_subscribers.begin(); it != m_subscribers.end(); ++it) {
            // check if the pointers agree
            if (it->first == id.data()) {
                m_subscribers.erase(it);
                return;
            }
        }
        assert_dbg(false, ExcUnknownSubscriberId(id, m_classname));
#endif
    }

    /** Get a subscription
     *
     * @param id The id to print if the subscription is not removed properly
     *           before deletion
     *
     * @note only has an effect if we are in DEBUG mode
     * */
    void subscribe(const std::string& id) const {
#ifdef DEBUG
        m_subscribers.push_front(std::make_pair(id.data(), id));

        // Set classname here, since this is actually executed by the
        // precise object we subscribe to and not the generic Subscribable
        // class. So here we have the "proper" type available in this.
        m_classname = std::string(typeid(*this).name());
#endif
    }

  private:
#ifdef DEBUG
    /** List to contain the pointer to the original string object, along
     *  with a copy of its actual content
     *
     * Marked as mutable in order to allow to subscribe / unsubscribe from
     * const references as well.
     */
    mutable std::list<std::pair<const char*, std::string>> m_subscribers;

    /**
     * Name of the actual child Subscribable class
     * Set on call of the subscribe function
     *
     * Marked as mutable in order to allow to subscribe / unsubscribe from
     * const references as well.
     * */
    mutable std::string m_classname{"(unknown)"};
#endif
};

}  // namespace linalgwrap

#endif  // SUBSCRIBABLE_H_
