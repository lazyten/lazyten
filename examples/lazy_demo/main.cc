//
// Copyright (C) 2017 by the lazyten authors
//
// This file is part of lazyten.
//
// lazyten is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// lazyten is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with lazyten. If not, see <http://www.gnu.org/licenses/>.
//

#include "timing.hh"
#include <iostream>
#include <lazyten/DiagonalMatrix.hh>
#include <lazyten/LazyMatrixWrapper.hh>
#include <lazyten/SmallMatrix.hh>
#include <lazyten/SmallVector.hh>
#include <lazyten/eigensystem.hh>
#include <lazyten/random.hh>
#include <lazyten/trans.hh>

using namespace lazyten;

// Define types
typedef double scalar_type;
typedef SmallMatrix<scalar_type> matrix_type;
typedef SmallVector<scalar_type> vector_type;
typedef DiagonalMatrix<matrix_type> diagmat_type;

int main() {
  //
  // Setup random matrices
  //
  const size_t size = 1000;
  diagmat_type diag(random<vector_type>(size));       // Diagonal matrix (lazy)
  matrix_type mat = random<matrix_type>(size, size);  // Stored matrix

  //
  // Build up an expression tree (no computation here)
  //
  auto start_tree = timing::now();

  auto sum = diag + mat;
  auto prod = trans(sum) * diag * sum;
  auto tree = mat + prod / 1e6;

  auto end_tree = timing::now();

  //
  // Apply to an vector
  //
  auto start_apply = timing::now();
  auto res = tree * random<vector_type>(size);
  auto end_apply = timing::now();

  //
  // Build a nice matrix
  //    (To add such that the diagonalisation becomes feasible)
  //
  vector_type diagonal(size);
  diagonal[0] = 11111;
  diagonal[1] = 54321;
  auto nice = make_diagmat(diagonal);

  //
  // Compute eigenpairs
  //
  auto start_sym = timing::now();
  auto tree_sym = trans(tree) + tree + nice;
  auto end_sym = timing::now();

  auto start_eigensys = timing::now();
  auto esolution = eigensystem_hermitian(tree_sym, 2, {{"which", "LR"}});
  auto end_eigensys = timing::now();

  //
  // Print results and timings
  //
  std::cout << "Result vector: " << res[0] << " " << res[1] << " " << res[2] << " ... "
            << res[size - 1] << std::endl;
  std::cout << "Result eigenvalues: " << esolution.evalues()[0] << " and "
            << esolution.evalues()[1] << std::endl
            << std::endl;

  std::cout << "Building the trees:   " << std::setw(6)
            << timing::time_in_ms(end_tree - start_tree) +
                     timing::time_in_ms(end_sym - start_sym)
            << " ms" << '\n'
            << "Applying the tree:    " << std::setw(6)
            << timing::time_in_ms(end_apply - start_apply) << " ms" << '\n'
            << "Solving eigensystem:  " << std::setw(6)
            << timing::time_in_ms(end_eigensys - start_eigensys) << " ms" << std::endl;

  return 0;
}
