#include <linalgwrap/DiagonalMatrix.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/eigensystem.hh>

using namespace linalgwrap;

template <typename Solution>
void print_solution(const Solution& solution) {
  for (size_t i = 0; i < solution.evalues().size(); ++i) {
    std::cout << "  (value: " << solution.evalues()[i]
              << "  ; vector: " << solution.evectors()[i] << ")" << std::endl;
  }
}

int main() {
  typedef double scalar_type;
  typedef SmallMatrix<scalar_type> matrix_type;

  // A hermitian test matrix
  matrix_type A{{76.52963291652406, -25.666510871673864, 34.43577632850656,
                 -45.36273708264996, 22.754426259474258, -11.270092801526317},
                {-25.666510871673864, -82.22395447280388, 51.60398579599948,
                 20.19904424373368, 46.73986164186124, -19.391628960288173},
                {34.43577632850656, 51.60398579599948, -44.14761987677002,
                 37.103640372970546, -23.036847741099365, 22.45800147301634},
                {-45.36273708264996, 20.19904424373368, 37.103640372970546,
                 52.40569142155442, 7.539858840849945, 28.047790646596525},
                {22.754426259474258, 46.73986164186124, -23.036847741099365,
                 7.539858840849945, -30.003084289104322, -80.44779473550247},
                {-11.270092801526317, -19.391628960288173, 22.45800147301634,
                 28.047790646596525, -80.44779473550247, 64.66122950159928}};

  // Compute the 4 largest real eigenpairs
  const size_t n_ep = 4;
  krims::ParameterMap map{{"which", "LR"}};

  //
  // Compute solution to the eigensystem using Arpack
  //
  map.update("method", "arpack");
  const auto solution_arpack = linalgwrap::eigensystem_hermitian(A, n_ep, map);
  std::cout << "Arpack eigenpairs: " << std::endl;
  print_solution(solution_arpack);
  std::cout << std::endl;

  //
  // Compute solution to the eigensystem using Arpack
  //
  map.update("method", "armadillo");
  const auto solution_armadillo = linalgwrap::eigensystem_hermitian(A, n_ep, map);

  std::cout << "Arpack eigenpairs: " << std::endl;
  print_solution(solution_armadillo);
  std::cout << std::endl;

  //
  // Compute solution for eigensystem given by
  // lazy matrix expression
  //
  map.update("method", "auto");
  SmallVector<scalar_type> diagonal{1, 2, -200, -1, 0, 0};
  auto diagmatrix = make_diagmat(std::move(diagonal));

  const auto sum = diagmatrix + A;
  const auto solution_lazy = linalgwrap::eigensystem_hermitian(sum, n_ep, map);
  std::cout << "Lazy matrix eigenpairs: " << std::endl;
  print_solution(solution_lazy);
  std::cout << std::endl;

  return 0;
}
