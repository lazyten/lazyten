#include "linalgwrap/ExceptionSystem.hh"

namespace linalgwrap {
// Set default effect of the assert_dbg macro to ABORT
ExceptionEffect AssertDbgEffect::m_eff = ExceptionEffect::ABORT;

void AssertDbgEffect::set(ExceptionEffect effect) { m_eff = effect; }
ExceptionEffect AssertDbgEffect::get() { return m_eff; }

}  // namespace linalgwrap
