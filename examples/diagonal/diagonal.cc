//
// Copyright (C) 2016 by the linalgwrap authors
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

#include "DiagonalUpdatable.hh"
#include <algorithm>
#include <iostream>
#include <linalgwrap/LazyMatrixWrapper.hh>
#include <linalgwrap/SmallMatrix.hh>
#include <linalgwrap/view.hh>

using namespace linalgwrap;

/** A genenaric rotation function, which works for all kind of
 *  stored or lazy rotation matrices or subject matrices
 *
 * Our dummy "algorithm" in this example.
 *  */
template <typename Matrix, typename RMatrix>
auto rotate(Matrix& m, RMatrix rotation_matrix)
      -> decltype(view::transpose(rotation_matrix) * m * rotation_matrix) {
    return view::transpose(rotation_matrix) * m * rotation_matrix;
}

int main() {
    // Define some types
    typedef double scalar_type;
    typedef SmallMatrix<scalar_type> stored_matrix_type;

    //
    // -------------------------
    //

    // Initialise a stored matrix:
    stored_matrix_type matrix{{42., 3.141592, 0.},   // First row
                              {-3.141592, -4., 0.},  // Second row
                              {0., 0., 1.}};         // Third row

    std::cout << "matrix = " << std::endl << matrix << std::endl << std::endl;

    // Initialise a diagonal matrix from a vector:
    SmallVector<scalar_type> diagonal{1.0, -4.0, 42.0};
    DiagonalUpdatable<stored_matrix_type> diag{diagonal};
    std::cout << "diag = " << std::endl << diag << std::endl << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //
    // Add matrices:
    //
    auto matrix_matrix = matrix + matrix;
    auto matrix_diag = matrix + diag;
    auto diag_diag = diag + diag;

    std::cout << "matrix+matrix = " << std::endl
              << matrix_matrix << std::endl
              << std::endl;

    std::cout << "matrix+diag = " << std::endl
              << matrix_diag << std::endl
              << std::endl;

    std::cout << "diag+diag = " << std::endl
              << diag_diag << std::endl
              << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //
    // Update the diagonal of the diag matrix
    //
    // This is done via a ParameterMap object, which maps an identifier string
    // to a pointer to an object. As pointer types both shared pointers as well
    // as SubscriptionPointers are allowed.
    // Here we use a SubscriptionPointer, which is implicitly constructed using
    // the make_subscription function from linalgwrap.
    ParameterMap map;
    diagonal = SmallVector<scalar_type>{-1., 2., 3.};
    map.update("diagonal", make_subscription(diagonal, "UpdateMap"));
    diag.update(map);

    // The output of the diag should change
    std::cout << "diag (updated diag) = " << std::endl
              << diag << std::endl
              << std::endl;

    // This will stay as it was.
    std::cout << "matrix+diag (updated diag) = " << std::endl
              << matrix_diag << std::endl
              << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    // Update the diag again, this time via matrix_diag
    diagonal = SmallVector<scalar_type>{4., 3., 2.};
    map.update("diagonal", make_subscription(diagonal, "UpdateMap"));
    matrix_diag.update(map);

    // This time this is unchanged
    std::cout << "diag (updated sum) = " << std::endl
              << diag << std::endl
              << std::endl;

    // Whereas this changes
    std::cout << "matrix+diag (updated sum) = " << std::endl
              << matrix_diag << std::endl
              << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //
    // Do a rotation, i.e. call the above "algorithm"
    //
    // Define a rotation matrix
    stored_matrix_type R{{0., 1., 0.},   // 1st row
                         {1., 0., 0.},   // 2nd row
                         {0., 0., 1.}};  // 3rd row

    std::cout << "The rotation matrix is" << std::endl
              << R << std::endl
              << std::endl;

    // Perform rotation:
    std::cout << "Rotation of matrix+diag = " << std::endl
              << rotate(matrix_diag, R) << std::endl
              << std::endl;
    std::cout << "Rotation of matrix = " << std::endl
              << rotate(matrix, R) << std::endl
              << std::endl;
    std::cout << "Rotation of diag = " << std::endl
              << rotate(diag, R) << std::endl
              << std::endl;

    std::cout << "------------------------------------------" << std::endl;

    //
    // Explicit cast to stored matrix
    //
    // Get a hard copy of the matrix+diag matrix expression
    // via performing an explicit cast to stored_matrix_type
    stored_matrix_type sum_copy = static_cast<stored_matrix_type>(matrix_diag);

    // Now do another update:
    diagonal = SmallVector<scalar_type>{-100., -100., -100.};
    map.update("diagonal", make_subscription(diagonal, "UpdateMap"));
    matrix_diag.update(map);

    // This of cause stays as it was, since it is now an evaluated
    // matrix in memory
    std::cout << "sum_copy (updated matrix+diag again) = " << std::endl
              << sum_copy << std::endl
              << std::endl;

    // This is the new, updated version.
    std::cout << "matrix+diag (updated matrix+diag again) = " << std::endl
              << matrix_diag << std::endl
              << std::endl;

    return 0;
}
