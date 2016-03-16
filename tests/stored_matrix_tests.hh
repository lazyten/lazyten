#pragma once
#include "matrix_tests.hh"

namespace linalgwrap {
namespace tests {
using namespace rc;

namespace stored_matrix_tests {

template <typename Matrix>
class TestingLibrary {
  public:
    typedef Matrix matrix_type;
    typedef typename matrix_type::size_type size_type;
    typedef typename matrix_type::scalar_type scalar_type;

    /** Run a sequence of rapidcheck checks in order to check
     * all basic functionality of a stored matrix */
    static bool run_check(std::string prefix = "");

  protected:
    // The testing library we use
    typedef matrix_tests::ComparativeTests<matrix_type, matrix_type> comptests;

    /** Helper-struct which effectively calls a test callable
     * from ComparativeTests using a randomly generated
     * matrix.
     */
    template <typename Callable>
    struct CallTest {
        Callable& call;
        void operator()() const;
        CallTest(Callable& c);
    };

    /** Helper-function to create a CallTest object */
    template <typename Callable>
    static CallTest<Callable> call_test(Callable& c);

    //
    // Extra tests for stored matrices
    //
    /** Test whether copying a stored matrix works */
    static void test_copy();

    /** Test whether setting elements via () works */
    static void test_setting_elements();

    /** Test whether setting elements via [] works */
    static void test_setting_elements_vectorised();

    /** Test whether setting elements via iterators works */
    static void test_read_write_iterator();
};

//
// ---------------------------------------
//

template <typename Matrix>
template <typename Callable>
void TestingLibrary<Matrix>::CallTest<Callable>::operator()() const {
    matrix_type m = *gen::arbitrary<matrix_type>().as("Matrix");
    call(m, m);
}

template <typename Matrix>
template <typename Callable>
TestingLibrary<Matrix>::CallTest<Callable>::CallTest(Callable& c)
      : call(c) {}

//
// TestingLibrary
//

template <typename Matrix>
template <typename Callable>
auto TestingLibrary<Matrix>::call_test(Callable& c) -> CallTest<Callable> {
    return CallTest<Callable>{c};
}

template <typename Matrix>
void TestingLibrary<Matrix>::test_copy() {
    matrix_type m = *gen::arbitrary<matrix_type>().as("Original matrix");
    matrix_type copy{m};

    // check that they are identical:
    RC_ASSERT(NumComp::is_equal_matrix(
          copy, m, std::numeric_limits<scalar_type>::epsilon()));
}

template <typename Matrix>
void TestingLibrary<Matrix>::test_setting_elements() {
    matrix_type m = *gen::arbitrary<matrix_type>().as("Matrix");

    auto modify_row =
          *gen::inRange<size_type>(0, m.n_rows()).as("Row to modify");
    auto modify_col =
          *gen::inRange<size_type>(0, m.n_cols()).as("Col to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    matrix_type m_copy{m};

    // Modify the value
    m(modify_row, modify_col) = value;

    // Check it:
    for (auto row : range(m.n_rows())) {
        for (auto col : range(m.n_cols())) {
            if (row == modify_row && col == modify_col) {
                RC_ASSERT(NumComp::is_equal(m(row, col), value));
            } else {
                RC_ASSERT(NumComp::is_equal(m(row, col), m_copy(row, col)));
            }
        }
    }
}

template <typename Matrix>
void TestingLibrary<Matrix>::test_setting_elements_vectorised() {
    matrix_type m = *gen::arbitrary<matrix_type>().as("Matrix");

    auto index = *gen::inRange<size_type>(0, m.n_rows() * m.n_cols())
                        .as("Index to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    // Keep a copy of the original:
    matrix_type m_copy{m};

    // Modify the value
    m[index] = value;

    // Check it:
    for (auto i : range(m.n_rows() * m.n_cols())) {
        if (i == index) {
            RC_ASSERT(NumComp::is_equal(m[i], value));
        } else {
            RC_ASSERT(NumComp::is_equal(m[i], m_copy[i]));
        }
    }
}

template <typename Matrix>
void TestingLibrary<Matrix>::test_read_write_iterator() {
    matrix_type m = *gen::arbitrary<matrix_type>().as("Matrix");

    auto modify_row =
          *gen::inRange<size_type>(0, m.n_rows()).as("Row to modify");
    auto modify_col =
          *gen::inRange<size_type>(0, m.n_cols()).as("Col to modify");
    auto value = *gen::arbitrary<scalar_type>().as("New value");

    for (auto it = m.begin(); it != m.end(); ++it) {
        if (it.row() == modify_row && it.col() == modify_col) {
            *it = value;
            break;
        }
    }

    auto it = std::begin(m);
    for (size_type i = 0; i < m.n_rows(); ++i) {
        for (size_type j = 0; j < m.n_cols(); ++j, ++it) {
            if (i == modify_row && j == modify_col) {
                RC_ASSERT(NumComp::is_equal(m(i, j), value));
            } else {
                RC_ASSERT(NumComp::is_equal(m(i, j), *it));
            }
        }
    }
}

template <typename Matrix>
bool TestingLibrary<Matrix>::run_check(std::string prefix) {
    bool res = true;

    // Test copying stored matrices
    res &= rc::check(prefix + "Test copying stored matrices", test_copy);

    // Read-only element access
    res &= rc::check(prefix + "Element access via () and []",
                     call_test(comptests::test_element_access));
    res &= rc::check(prefix + "Element access via extract_block",
                     call_test(comptests::test_extract_block));
    res &= rc::check(prefix + "Data access via add_block_to",
                     call_test(comptests::test_add_block_to));
    res &= rc::check(prefix + "Read-only iterator of small matrices",
                     call_test(comptests::test_iterator));

    // Read-write element access
    res &=
          rc::check(prefix + "Altering elements via ()", test_setting_elements);
    res &= rc::check(prefix + "Altering elements via []",
                     test_setting_elements_vectorised);
    res &= rc::check(prefix + "Altering elements via iterator",
                     test_read_write_iterator);

    // Operations
    res &= rc::check(prefix + "Multiplication by scalar",
                     call_test(comptests::test_mutiply_scalar));
    res &= rc::check(prefix + "Divide by scalar",
                     call_test(comptests::test_divide_scalar));
    res &= rc::check(prefix + "Add a matrix",
                     call_test(comptests::template test_add<Matrix>));
    res &= rc::check(prefix + "Subtract a matrix",
                     call_test(comptests::template test_subtract<Matrix>));
    res &= rc::check(prefix + "Matrix multiplication",
                     call_test(comptests::template test_multiply_by<Matrix>));

    return res;
}

}  // namespace stored_matrix_tests
}  // namespace tests
}  // namespace linalgwrap
