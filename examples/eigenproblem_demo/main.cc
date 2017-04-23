//
// Copyright (C) 2017 by the linalgwrap authors
//
// This file is part of linalgwrap.
//
// linalgwrap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// linalgwrap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with linalgwrap. If not, see <http://www.gnu.org/licenses/>.
//

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
  krims::GenMap map{{"which", "LR"}};

  //
  // Compute solution to the eigensystem using Arpack
  //
  map.update("method", "arpack");
  const auto solution_arpack = linalgwrap::eigensystem_hermitian(A, n_ep, map);
  std::cout << "Arpack eigenpairs: " << std::endl;
  print_solution(solution_arpack);
  std::cout << std::endl;

  //
  // Compute solution to the eigensystem using Armadillo
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
  {
    map.update("method", "auto");
    SmallVector<scalar_type> diagonal{1, 2, -200, -1, 0, 0};
    auto diagmatrix = make_diagmat(std::move(diagonal));

    const auto sum = diagmatrix + A;
    const auto solution_lazy = linalgwrap::eigensystem_hermitian(sum, n_ep, map);
    std::cout << "Lazy matrix eigenpairs: " << std::endl;
    print_solution(solution_lazy);
    std::cout << std::endl;
  }

  //
  // Compute solution for generalised eigensystem
  // given by two lazy matrix expressions
  //
  {
    map.update("method", "auto");
    SmallMatrix<scalar_type> A{
          {-37.827087025623456, 7.737383999517007, -6.868688554610827, 0.4582402410539501,
           20.25022267287818, -7.264272891184383, 17.13528963873486, -0.6958437952983445},
          {7.737383999517007, -8.907240678731085, -10.491465613577162,
           -32.400293305571346, 1.9050289475850448, 38.08086976760322, 33.58637188467661,
           16.93991633946149},
          {-6.868688554610827, -10.491465613577162, -40.64686080355078,
           26.642223861646016, 26.870722616736714, 10.117350139698132,
           -0.19698340817648585, 14.714415823222282},
          {0.4582402410539501, -32.400293305571346, 26.642223861646016,
           -28.78493388224115, -10.56060566681876, -7.7461795092692824,
           -9.341941942209004, -3.682243002650196},
          {20.25022267287818, 1.9050289475850448, 26.870722616736714, -10.56060566681876,
           -35.21142298086606, 4.5353204853178966, 41.307563515495005,
           -10.243149796300195},
          {-7.264272891184383, 38.08086976760322, 10.117350139698132, -7.7461795092692824,
           4.5353204853178966, 46.221528657079915, -11.029394211181724,
           33.71103101586654},
          {17.13528963873486, 33.58637188467661, -0.19698340817648585, -9.341941942209004,
           41.307563515495005, -11.029394211181724, 24.566705894503286,
           -5.017034388900839},
          {-0.6958437952983445, 16.93991633946149, 14.714415823222282, -3.682243002650196,
           -10.243149796300195, 33.71103101586654, -5.017034388900839,
           13.035884650583398}};

    SmallVector<scalar_type> diagonal{1, 2, 3, 1, 2, 4, 8, 9};
    auto diagmatrix = make_diagmat(std::move(diagonal));

    const size_t n_ep = 2;
    const auto sum = A + 100 * diagmatrix;
    const auto solution_lgen =
          linalgwrap::eigensystem_hermitian(sum, diagmatrix, n_ep, map);
    std::cout << "Lazy matrix generalised eigenpairs: " << std::endl;
    print_solution(solution_lgen);
    std::cout << std::endl;
  }

  return 0;
}
