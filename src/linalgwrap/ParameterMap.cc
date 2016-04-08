#include "linalgwrap/ParameterMap.hh"

namespace linalgwrap {

//
// Entry subclass
//
ParameterMap::Entry::Entry()
      : m_object_ptr(nullptr), m_via_subscription_ptr(false) {
#ifdef DEBUG
    m_type_name = "";
#endif
}

//
// ParameterMap
//

/** Remove an element */
void ParameterMap::erase(const std::string& key) {
    assert_dbg(exists(key), ExcUnknownKey(key));
    m_container.erase(key);
}

/** Check weather a key exists */
bool ParameterMap::exists(const std::string& key) const {
    return m_container.find(key) != std::end(m_container);
}

}  // linalgwrap
