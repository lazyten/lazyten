#pragma once
#include "matrix_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

/** Namespace for default tests for stored matrices */
namespace stored_matrix_tests {

template <typename Matrix>
class TestingLibrary {
  public:
    typedef Matrix matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    /** Construct a testing library instance in order to check
     *  all basic functionality of a stored matrix
     *
     *  \param prefix A prefix to use in the rapidcheck description string
     *  \param tolerance Numeric tolerance used for comparison in (some) checks.
     *                   The checks for which machine epsilon is used as
     *                   tolerance are:
     *  - Copying a stored matrix
     *  - Setting and getting elements via [], () or iterator
     */
    TestingLibrary(std::string prefix = "",
                   double tolerance = TestConstants::default_num_tol);

    void run_checks() const;

  private:
    // The testing library we use
    typedef matrix_tests::ComparativeTests<matrix_type, matrix_type> comptests;

    // The caller type we use
    typedef typename comptests::template RapidcheckTestableGenerator<
          matrix_type> callgen_type;

    /** Identity operation on the matrix */
    static constexpr matrix_type identity(matrix_type m) { return m; };

    /** Argument generation */
    static constexpr matrix_type argsgen() {
        return *gen::arbitrary<matrix_type>().as("Matrix");
    };

    std::string m_prefix;
    callgen_type m_gen;
};

//
// ---------------------------------------
//

template <typename Matrix>
TestingLibrary<Matrix>::TestingLibrary(std::string prefix, double tolerance)
      : m_prefix{prefix}, m_gen{argsgen, identity, identity, tolerance} {}

template <typename Matrix>
void TestingLibrary<Matrix>::run_checks() const {
    const double eps = std::numeric_limits<scalar_type>::epsilon();

    // Test copying stored matrices
    CHECK(rc::check(m_prefix + "Test copying stored matrices",
                    m_gen.generate(comptests::test_copy, eps)));

    // Read-only element access
    CHECK(rc::check(m_prefix + "Element access via () and []",
                    m_gen.generate(comptests::test_element_access, eps)));
    CHECK(rc::check(m_prefix + "Element access via extract_block",
                    m_gen.generate(comptests::test_extract_block)));
    CHECK(rc::check(m_prefix + "Data access via add_block_to",
                    m_gen.generate(comptests::test_add_block_to)));
    CHECK(rc::check(m_prefix + "Read-only iterator of small matrices",
                    m_gen.generate(comptests::test_readonly_iterator, eps)));

    // Read-write element access
    CHECK(rc::check(
          m_prefix + "Altering elements via ()",
          m_gen.generate(comptests::test_setting_elements_indexed, eps)));
    CHECK(rc::check(
          m_prefix + "Altering elements via []",
          m_gen.generate(comptests::test_setting_elements_vectorised, eps)));
    CHECK(rc::check(m_prefix + "Altering elements via iterator",
                    m_gen.generate(comptests::test_readwrite_iterator, eps)));

    // Operations
    CHECK(rc::check(m_prefix + "Multiplication by scalar",
                    m_gen.generate(comptests::test_mutiply_scalar)));
    CHECK(rc::check(m_prefix + "Divide by scalar",
                    m_gen.generate(comptests::test_divide_scalar)));
    CHECK(rc::check(m_prefix + "Add a matrix",
                    m_gen.generate(comptests::template test_add<Matrix>)));
    CHECK(rc::check(m_prefix + "Subtract a matrix",
                    m_gen.generate(comptests::template test_subtract<Matrix>)));
    CHECK(rc::check(
          m_prefix + "Matrix multiplication",
          m_gen.generate(comptests::template test_multiply_by<Matrix>)));
}

}  // namespace stored_matrix_tests
}  // namespace tests
}  // namespace linalgwrap
