#ifndef LINALGWRAP_PARAMETER_WRAP_HH_
#define LINALGWRAP_PARAMETER_WRAP_HH_

#pragma once
#include "linalgwrap/Exceptions.hh"
#include "linalgwrap/SubscriptionPointer.hh"
#include <map>
#include <memory>
#include <string>

namespace linalgwrap {

/** ParameterMap is essentially a map from std::string to
 *  a shared pointer of void. This way an arbitrary amount of arbitrary
 *  objects can be passed around using this object and they can be quickly
 *  accessed by the std::string key.
 */
class ParameterMap {
  public:
    //
    // Exception declaration:
    //
    /** Exception to indicate that a wrong type was requested */
    DefException3(ExcWrongTypeRequested, std::string, std::string, std::string,
                  << "Requested type " << arg1
                  << " from ParameterMap using key " << arg2
                  << ". The value however has type " << arg3 << ".");

    /** Exception thrown if a key is not valid */
    DefException1(ExcUnknownKey, std::string, << "The key " << arg1
                                              << " is unknown.");

  private:
    //
    // The ExtractValueFromPointer class
    //
    /** Helper class to make sure that no compiler error gets produced
     * when T is not a Subscribable
     */
    template <typename T,
              typename = typename std::is_base_of<Subscribable, T>::type>
    struct DereferencePointerSubscriptionPointer {
        typedef T& result_type;
        typedef const std::shared_ptr<void> argument_type;

        T& operator()(std::shared_ptr<void> p) const;
    };

    /** Specialisation for objects which are not subscribable
     *
     * Does nothing.
     * */
    template <typename T>
    struct DereferencePointerSubscriptionPointer<T, std::false_type> {
        typedef T& result_type;
        typedef const std::shared_ptr<void> argument_type;

        T& operator()(std::shared_ptr<void>);

      private:
        /** A dummy m_t to return */
        T m_t;
    };

    //
    // The entry class
    //
    class Entry {
      public:
        /** Default constructor: Construct empty object */
        Entry();

        template <typename T>
        explicit Entry(std::shared_ptr<T> ptr);

        template <typename T>
        explicit Entry(SubscriptionPointer<T> ptr);

        template <typename T>
        T& get(const std::string& key);

        template <typename T>
        const T& get(const std::string& key) const;

      private:
        // The stored pointer
        std::shared_ptr<void> m_object_ptr;

        // Is the object_ptr a pointer to the object
        // or a pointer to a subscription_pointer
        bool m_via_subscription_ptr;
#ifdef DEBUG
        std::string m_type_name;
#endif
    };

  public:
    /** Insert or update using a shared pointer */
    template <typename T>
    void update(std::string key, std::shared_ptr<T> object_ptr);

    /** Insert or update using a SubscriptionPointer */
    template <typename T>
    void update(std::string key, SubscriptionPointer<T> object_ptr);

    /** Insert or update a key with a copy of an element */
    template <typename T>
    void update_copy(std::string key, T object);

    /** Remove an element */
    void erase(const std::string& key);

    /** Check weather a key exists */
    bool exists(const std::string& key) const;

    /** Get the value of an element
     */
    template <typename T>
    T& at(const std::string& key);

    /** Return the value at a given key in a specific type
     */
    template <typename T>
    const T& at(const std::string& key) const;

  private:
    std::map<std::string, Entry> m_container;
};

//
// -----------------------------------------------------------------
//

//
// FullDereference
//
template <typename T, typename B>
T& ParameterMap::DereferencePointerSubscriptionPointer<T, B>::operator()(
      std::shared_ptr<void> p) const {
    // Extract pointer to subscription pointer
    auto ptr_ptr = std::static_pointer_cast<SubscriptionPointer<T>>(p);

    // return the value by twice dereferencing:
    return *(*ptr_ptr);
}

template <typename T>
T& ParameterMap::DereferencePointerSubscriptionPointer<T, std::false_type>::
operator()(std::shared_ptr<void>) {
    assert_dbg(false, ExcNotImplemented());
    return m_t;
}

//
// Entry subclass
//
template <typename T>
ParameterMap::Entry::Entry(std::shared_ptr<T> ptr)
      : m_object_ptr(ptr), m_via_subscription_ptr(false) {
#ifdef DEBUG
    m_type_name = typeid(T).name();
#endif
}

template <typename T>
ParameterMap::Entry::Entry(SubscriptionPointer<T> ptr)
      : m_object_ptr(std::make_shared<SubscriptionPointer<T>>(ptr)),
        m_via_subscription_ptr(true) {
#ifdef DEBUG
    m_type_name = typeid(T).name();
#endif
}

template <typename T>
T& ParameterMap::Entry::get(const std::string& key) {
    // check that the correct type is requested:
    assert_dbg(m_type_name == typeid(T).name(),
               ExcWrongTypeRequested(typeid(T).name(), key, m_type_name));

    if (m_via_subscription_ptr) {
        DereferencePointerSubscriptionPointer<T> dereference;
        return dereference(m_object_ptr);
    } else {
        // Extract pointer to type and dereference it
        return *std::static_pointer_cast<T>(m_object_ptr);
    }
}

template <typename T>
const T& ParameterMap::Entry::get(const std::string& key) const {
    // check that the correct type is requested:
    assert_dbg(m_type_name == typeid(T).name(),
               ExcWrongTypeRequested(typeid(T).name(), key, m_type_name));

    if (m_via_subscription_ptr) {
        DereferencePointerSubscriptionPointer<T> dereference;
        return dereference(m_object_ptr);
    } else {
        // Extract pointer to type and dereference it
        return *std::static_pointer_cast<T>(m_object_ptr);
    }
}

//
// ParameterMap
//
template <typename T>
void ParameterMap::update(std::string key, std::shared_ptr<T> object_ptr) {
    // Insert or update a new element:
    m_container[key] = Entry{object_ptr};
}

/** Insert or update using a SubscriptionPointer */
template <typename T>
void ParameterMap::update(std::string key, SubscriptionPointer<T> object_ptr) {
    // Insert or update a new element:
    m_container[key] = Entry{object_ptr};
}

/** Insert or update a key with a copy of an element */
template <typename T>
void ParameterMap::update_copy(std::string key, T object) {
    m_container[key] = Entry{std::make_shared<T>(object)};
}

/** Get the value of an element
 */
template <typename T>
T& ParameterMap::at(const std::string& key) {
    assert_dbg(exists(key), ExcUnknownKey(key));
    return m_container.at(key).get<T>(key);
}

/** Return the value at a given key in a specific type
 */
template <typename T>
const T& ParameterMap::at(const std::string& key) const {
    assert_dbg(exists(key), ExcUnknownKey(key));
    return m_container.at(key).get<T>(key);
}

}  // linalgwrap
#endif
