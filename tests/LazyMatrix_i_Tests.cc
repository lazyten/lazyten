#include "lazy_matrix_tests.hh"
#include <catch.hpp>
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/LazyMatrix_i.hh>
#include <rapidcheck.h>

namespace linalgwrap {
namespace tests {
using namespace rc;

namespace lazymatrix_i_tests {

template <typename StoredMatrix>
class SimpleLazyMatrix : public LazyMatrix_i<StoredMatrix> {
  public:
    typedef LazyMatrix_i<StoredMatrix> base_type;
    typedef typename base_type::stored_matrix_type stored_matrix_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::scalar_type scalar_type;
    typedef typename base_type::lazy_matrix_expression_ptr_type
          lazy_matrix_expression_ptr_type;

    SimpleLazyMatrix(stored_matrix_type m) : m_stored{m} {}

    size_type n_rows() const override { return m_stored.n_rows(); }
    size_type n_cols() const override { return m_stored.n_cols(); }
    scalar_type operator()(size_type i, size_type j) const override {
        return m_stored(i, j);
    }

    lazy_matrix_expression_ptr_type clone() const override {
        return lazy_matrix_expression_ptr_type(new SimpleLazyMatrix(*this));
    }

  private:
    stored_matrix_type m_stored;
};
}

TEST_CASE("LazyMatrix_i abstract class", "[LazyMatrix_i]") {
    using namespace lazymatrix_i_tests;

    // Make sure that the program does not get aborted
    exceptions::assert_dbg_effect = exceptions::ExceptionEffect::THROW;

    // Test constructor
    // Test swapping

    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;
    typedef SimpleLazyMatrix<stored_matrix_type> lazy_matrix_type;
    typedef typename stored_matrix_type::size_type size_type;

    // Generator for the args
    auto args_generator = [] {
        stored_matrix_type m =
              *gen::arbitrary<stored_matrix_type>().as("Inner matrix");
        // TODO remove to test empty matrix case.
        RC_PRE(m.n_cols() > size_type{0} && m.n_rows() > size_type{0});
        return m;
    };

    // Generator for the model.
    auto model_generator = [](stored_matrix_type m) { return m; };

    // Generator for the sut
    auto lazy_generator = [](stored_matrix_type m) {
        lazy_matrix_type wrap{m};
        return wrap;
    };

    SECTION("Default lazy matrix tests") {
        typedef lazy_matrix_tests::TestingLibrary<lazy_matrix_type,
                                                  decltype(args_generator())>
              testlib;

        testlib lib{args_generator, lazy_generator, model_generator,
                    "LazyMatrix_i: ", TestConstants::default_num_tol};
        lib.run_checks();
    }
}

}  // tests
}  // linalgwrap
