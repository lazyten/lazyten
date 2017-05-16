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

#include "matrices.hh"
#include <linalgwrap/DiagonalMatrix.hh>
#include <linalgwrap/SmallVector.hh>
#include <linalgwrap/eigensystem.hh>

template <typename Solution>
void print_solution(const Solution& solution) {
  for (size_t i = 0; i < solution.evalues().size(); ++i) {
    std::cout << "  (value: " << solution.evalues()[i]
              << "  ; vector: " << solution.evectors()[i] << ")" << std::endl;
  }
}

int main() {
  using namespace linalgwrap;

  // Compute the 4 largest real eigenpairs of the test matrix a
  const size_t n_ep = 4;
  krims::GenMap map{{"which", "LR"}};

  //
  // Compute solution to the eigensystem using Arpack
  //
  map.update("method", "arpack");
  const auto solution_arpack = linalgwrap::eigensystem_hermitian(mat_a, n_ep, map);
  std::cout << "Arpack eigenpairs: " << std::endl;
  print_solution(solution_arpack);
  std::cout << std::endl;

  //
  // Compute solution to the eigensystem using Lapack
  //
  map.update("method", "lapack");
  const auto solution_armadillo = linalgwrap::eigensystem_hermitian(mat_a, n_ep, map);

  std::cout << "Lapack eigenpairs: " << std::endl;
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

    const auto sum = diagmatrix + mat_a;
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
    SmallVector<scalar_type> diagonal{1, 2, 3, 1, 2, 4, 8, 9};
    auto diagmatrix = make_diagmat(std::move(diagonal));
    const auto sum = mat_b + 100 * diagmatrix;

    // Compute 2 largest eigenpairs
    const size_t n_ep = 2;
    const krims::GenMap params{{"method", "auto"}, {"which", "LR"}};
    const auto solution_lgen =
          linalgwrap::eigensystem_hermitian(sum, diagmatrix, n_ep, params);

    std::cout << "Lazy matrix generalised eigenpairs: " << std::endl;
    print_solution(solution_lgen);
    std::cout << std::endl;
  }

  return 0;
}
