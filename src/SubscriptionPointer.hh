#ifndef SUBSCRIPTIONPOINTER_H_
#define SUBSCRIPTIONPOINTER_H_

#include <utility>
#include <string>
#include <memory>
#include <type_traits>
#include "Subscribable.hh"

namespace linalgwrap {

// Forward-declare Subscribable
class Subscribable;

template <typename T>
class SubscriptionPointer {
    static_assert(std::is_base_of<Subscribable, T>::value,
                  "T must be a child class of Subscribable");

  public:
    /** A swap function for Subscription pointers */
    friend void swap(SubscriptionPointer& first, SubscriptionPointer& second) {
        using std::swap;

        swap(first.m_subscriber_id_ptr, second.m_subscriber_id_ptr);
        swap(first.m_subscribed_obj_ptr, second.m_subscribed_obj_ptr);
    }

    /** \brief Create an empty Subscription pointer pointing to no object at all
     *
     * Actual subscription to an object may be done using the reset function.
     *
     * \param subscriber_id Id used for the subscription.
     * */
    explicit SubscriptionPointer(const std::string& subscriber_id)
          : m_subscriber_id_ptr(new std::string(subscriber_id)),
            m_subscribed_obj_ptr(nullptr) {}

    /** Create a Subscription pointer from a Subscribable Class and an
     * identifier
     *
     * \param subscriber_id  identifier used for subscription
     * \param s object to subscribe to
     * */
    SubscriptionPointer(const std::string& subscriber_id, T& s)
          : SubscriptionPointer(subscriber_id) {
        // Register this subscription
        register_at(&s);
    }

    /** Destructor */
    ~SubscriptionPointer() {
        // Unregister from everything
        register_at(nullptr);
    }

    /** Copy constructor */
    SubscriptionPointer(const SubscriptionPointer& other)
          : m_subscriber_id_ptr(new std::string(*(other.m_subscriber_id_ptr))),
            m_subscribed_obj_ptr(nullptr) {

        // Register the subscription
        register_at(other.m_subscribed_obj_ptr);
    }

    /** Move constructor */
    SubscriptionPointer(SubscriptionPointer&& other)
          : m_subscriber_id_ptr(new std::string(*other.m_subscriber_id_ptr)),
            m_subscribed_obj_ptr(nullptr) {
        // Register us:
        register_at(other.m_subscribed_obj_ptr);

        // Unregister the other object:
        other.register_at(nullptr);
    }

    /** Assignment operator */
    SubscriptionPointer& operator=(SubscriptionPointer other) {
        swap(*this, other);
        return *this;
    }

    /** \brief Reset the subscription and subscribe to a new object.
     *
     * Unsubscribes from the old object and subscribes to the new object
     * with the exact same subscriber id
     *
     * \param s The new object to subscribe to
     * */
    void reset(T& s) { register_at(&s); }

    /** \brief Reset the subscription and become empty */
    void reset() { register_at(nullptr); }

    /** \brief Check if this object is empty or not */
    explicit operator bool() const { return m_subscribed_obj_ptr != nullptr; }

    /** \brief Raw access to the inner pointer */
    T* get() const { return m_subscribed_obj_ptr; }

    /** Dereference object */
    T& operator*() const { return *m_subscribed_obj_ptr; }

    /** Dereference object member */
    T* operator->() const { return m_subscribed_obj_ptr; }

    const std::string& subscriber_id() { return *m_subscriber_id_ptr; }

  private:
    /** Register at the given object */
    void register_at(T* object_ptr) {
        if (m_subscribed_obj_ptr) {
            // Unregister subscription of this pointer from old object
            m_subscribed_obj_ptr->unsubscribe(m_subscriber_id_ptr);
            m_subscribed_obj_ptr = nullptr;
        }

        if (object_ptr) {
            // Register subscription at new object
            object_ptr->subscribe(m_subscriber_id_ptr);
            m_subscribed_obj_ptr = object_ptr;
        }
    }

    // Use a unique_ptr of string here.
    std::shared_ptr<const std::string> m_subscriber_id_ptr;
    T* m_subscribed_obj_ptr;
};

// TODO improve subscribable - subscription system and remove
//      some of the problems mentioned below.
/**
 * \brief Convenience wrapper for making subscription pointers
 * and subscribing to objects
 *
 * When making subscriptions to a reference in this
 * way, the mechanism protects you from deallocating
 * the referenced object and leaving a dangling pointer
 * behind.
 *
 * It does *not* protect you from doing silly things to
 * the actual memory, however. So if you move the object
 * somewhere else, i.e. invalidate its memory while
 * still holding subscriptions to it, this mechanism
 * does not detect and complain about this.
 *
 * Similarly if you make changes to the object after
 * having made these subscriptions, these changes are
 * of course also seen from the Subscription pointers
 * which hold a subscription to said object.
 */
template <typename T>
inline SubscriptionPointer<T> make_subscription(
      T& object, const std::string& subscriber_id) {
    return SubscriptionPointer<T>(subscriber_id, object);
}

//
// == and != comparison operators:
//
/** Compare if the SubscriptionPointers point to the same object */
template <typename T>
inline bool operator==(const SubscriptionPointer<T>& lhs,
                       const SubscriptionPointer<T>& rhs) {
    return lhs.get() == rhs.get();
}

/** Compare if the SubscriptionPointers do not point to the same object */
template <typename T>
inline bool operator!=(const SubscriptionPointer<T>& lhs,
                       const SubscriptionPointer<T>& rhs) {
    return !operator==(lhs, rhs);
}

/** Compare if the Subscription pointer points to no object, i.e. stores a
 * nullptr */
template <typename T>
inline bool operator==(const SubscriptionPointer<T>& lhs, std::nullptr_t) {
    return lhs.get() == nullptr;
}

/** Compare if the Subscription pointer points to no object, i.e. stores a
 * nullptr */
template <typename T>
inline bool operator==(std::nullptr_t, const SubscriptionPointer<T>& rhs) {
    return rhs == nullptr;
}

/** Compare if the Subscription pointer does not point to no object, i.e. stores
 * no nullptr */
template <typename T>
inline bool operator!=(const SubscriptionPointer<T>& lhs, std::nullptr_t) {
    return !operator==(lhs, nullptr);
}

/** Compare if the Subscription pointer does not point to no object, i.e. stores
 * no nullptr */
template <typename T>
inline bool operator!=(std::nullptr_t, const SubscriptionPointer<T>& rhs) {
    return !operator==(nullptr, rhs);
}

}  // namespace linalgwrap

#endif  // SUBSCRIPTIONPOINTER_H_
