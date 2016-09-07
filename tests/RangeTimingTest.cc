#include <chrono>
#include <linalgwrap/Range.hh>
#include <vector>

static const size_t size = 10 * 1000 * 1000;

int main() {
    double a = 2.3;
    double b = 4.5;
    double res;

    auto p1 = std::chrono::system_clock::now();
    for (size_t i = 0; i < size; ++i) {
        res = a * b * i;
    }
    auto p2 = std::chrono::system_clock::now();
    for (size_t i : linalgwrap::range(size)) {
        res = a * b * i;
    }
    auto p3 = std::chrono::system_clock::now();

    auto duration_normal = p2 - p1;
    auto duration_range = p3 - p2;

    std::cout << "normal: "
              << duration_normal.count() / static_cast<double>(size)
              << std::endl;
    std::cout << "range:  "
              << duration_range.count() / static_cast<double>(size)
              << std::endl;
    return 0;
}
