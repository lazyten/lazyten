//
// Copyright (C) 2016-17 by the linalgwrap authors
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

#pragma once
#include "BlockDiagonalMatrix.hh"
#include "EigensystemSolver.hh"
#include <iterator>

namespace linalgwrap {

/** \brief Solve a normal hermitian eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 **/
template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value &&
                                    !IsBlockDiagonalMatrix<Matrix>::value,
                              size_t>::type n_ep = Constants<size_t>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a Hermitian eigensystem in blocked form
 *
 * The problem is essentially split up in a problem over the
 * individual blocks. The obtained eigenvectors are padded
 * with explicit zeros, such that they are only acting the
 * block from which they were obtained. The order is exactly
 * the order obtained from the individual eigensolvers,
 * pasted together block by block.
 *
 * \param n_ep Total number of eigenpairs to obtain.
 *             The number will be evenly spread between
 *             the blocks with earlier blocks getting a larger
 *             part of the share to break the tie.
 *
 * \note The interface of this function is very likely to change
 *       again in the near future. Especially the meaning
 *       of n_ep could change.
 */
template <typename BlockMatrix>
EigensolutionTypeFor<true, BlockMatrix> eigensystem_hermitian(
      const BlockMatrix& A,
      typename std::enable_if<IsBlockDiagonalMatrix<BlockMatrix>::value, size_t>::type
            n_ep = Constants<size_t>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a generalised hermitian eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value &&
                                    !IsBlockDiagonalMatrix<MatrixA>::value &&
                                    !IsBlockDiagonalMatrix<MatrixB>::value,
                              size_t>::type n_ep = Constants<size_t>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a generalised Hermitian eigensystem in blocked form
 *
 * The problem is essentially split up in a problem over the
 * individual blocks. The obtained eigenvectors are padded
 * with explicit zeros, such that they are only acting the
 * block from which they were obtained. The order is exactly
 * the order obtained from the individual eigensolvers,
 * pasted together block by block.
 *
 * \param n_ep Total number of eigenpairs to obtain.
 *             The number will be evenly spread between
 *             the blocks with earlier blocks getting a larger
 *             part of the share to break the tie.
 *
 * \note The interface of this function is very likely to change
 *       again in the near future. Especially the meaning
 *       of n_ep could change.
 */
template <typename BlockMatrixA, typename BlockMatrixB>
EigensolutionTypeFor<true, BlockMatrixA> eigensystem_hermitian(
      const BlockMatrixA& A, const BlockMatrixB& B,
      typename std::enable_if<IsBlockDiagonalMatrix<BlockMatrixA>::value &&
                                    IsBlockDiagonalMatrix<BlockMatrixB>::value,
                              size_t>::type n_ep = Constants<size_t>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a normal eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. The number of understood options and
 * their default values and supported values depend on the underlying
 * eigensolver which is used. The options which are supported for all
 * eigensolvers can be found in the documentation of EigensystemSolver and
 * in the struct EigensystemSolverKeys.
 *
 * \param A     Matrix to diagonalise
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep = Constants<typename Matrix::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

/** Solve a generalised eigensystem
 *
 * The Parameter map \t map may be used to provide configurable parameters
 * to the underlying Eigensolver. See the non-generalised version of
 * eigensolver_hermitian for which parameters are usually supported
 * and there default values are.
 *
 * \param A     Matrix to diagonalise
 * \param B     Metric to use
 * \param n_ep  Number of eigenpairs to obtain
 *              (Constants<typename Matrix::size_type>::all for all)
 * \param map   Specify some solver parameters
 *
 * \throws      Subclass of SolverException in case there is an error.
 **/
template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep =
            Constants<typename MatrixA::size_type>::all,
      const krims::GenMap& map = krims::GenMap());

//
// ------------------------------------------------------
//

template <typename Matrix>
EigensolutionTypeFor<true, Matrix> eigensystem_hermitian(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value &&
                                    !IsBlockDiagonalMatrix<Matrix>::value,
                              size_t>::type n_ep,
      const krims::GenMap& map) {
  typedef Eigenproblem<true, Matrix> problem_type;
  problem_type problem{A, n_ep};

  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename BlockMatrix>
EigensolutionTypeFor<true, BlockMatrix> eigensystem_hermitian(
      const BlockMatrix& A,
      typename std::enable_if<IsBlockDiagonalMatrix<BlockMatrix>::value, size_t>::type
            n_ep,
      const krims::GenMap& map) {
  typedef typename EigensolutionTypeFor<true, BlockMatrix>::evector_type evector_type;
  //  typedef typename EigensolutionTypeFor<true, BlockMatrix>::evalue_type evalue_type;

  // How many eigenpairs to compute per block
  const size_t n_ep_per_block = n_ep / A.diag_blocks().size();
  const size_t n_ep_rest = n_ep - n_ep_per_block * A.diag_blocks().size();

  EigensolutionTypeFor<true, BlockMatrix> ret;

  // Solve individual problems:
  size_t begin_index = 0;  // of the current block
  for (auto it = std::begin(A.diag_blocks()); it != std::end(A.diag_blocks()); ++it) {
    // How many eigenpairs should we solve for in this block?
    const size_t b_idx = it - std::begin(A.diag_blocks());
    const size_t n_ep_block = n_ep_per_block + (b_idx < n_ep_rest ? 1 : 0);

    auto block_soln = eigensystem_hermitian(*it, n_ep_block, map);
    assert_dbg(block_soln.n_ep() == n_ep_block, krims::ExcInternalError());

    // Lambda to pad the vector with zeros at beginning and end
    auto pad_with_zeros = [&begin_index, &A](const evector_type& v) {
      evector_type padded(A.n_cols());
      std::copy(v.begin(), v.end(), padded.begin() + begin_index);
      return padded;
    };

    // Copy into solution data structure and pad the eigenvectors with zeros
    std::copy(block_soln.evalues().begin(), block_soln.evalues().end(),
              std::back_inserter(ret.evalues()));
    std::transform(block_soln.evectors().begin(), block_soln.evectors().end(),
                   std::back_inserter(ret.evectors()), pad_with_zeros);

    begin_index += it->n_rows();
  }
  assert_dbg(begin_index == A.n_rows(), krims::ExcInternalError());

  return ret;
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<true, MatrixA> eigensystem_hermitian(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value &&
                                    !IsBlockDiagonalMatrix<MatrixA>::value &&
                                    !IsBlockDiagonalMatrix<MatrixB>::value,
                              size_t>::type n_ep,
      const krims::GenMap& map) {
  typedef Eigenproblem<true, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename BlockMatrixA, typename BlockMatrixB>
EigensolutionTypeFor<true, BlockMatrixA> eigensystem_hermitian(
      const BlockMatrixA& A, const BlockMatrixB& B,
      typename std::enable_if<IsBlockDiagonalMatrix<BlockMatrixA>::value &&
                                    IsBlockDiagonalMatrix<BlockMatrixB>::value,
                              size_t>::type n_ep,
      const krims::GenMap& map) {
  assert_size(A.diag_blocks().size(), B.diag_blocks().size());

  typedef typename EigensolutionTypeFor<true, BlockMatrixA>::evector_type evector_type;
  //  typedef typename EigensolutionTypeFor<true, BlockMatrixA>::evalue_type evalue_type;

  // How many eigenpairs to compute per block
  const size_t n_ep_per_block = n_ep / A.diag_blocks().size();
  const size_t n_ep_rest = n_ep - n_ep_per_block * A.diag_blocks().size();

  EigensolutionTypeFor<true, BlockMatrixA> ret;

  // Solve individual problems:
  size_t begin_index = 0;  // of the current block
  auto itb = std::begin(B.diag_blocks());
  for (auto ita = std::begin(A.diag_blocks()); ita != std::end(A.diag_blocks());
       ++ita, ++itb) {
    assert_size(ita->size(), itb->size());

    // How many eigenpairs should we solve for in this block?
    const size_t b_idx = ita - std::begin(A.diag_blocks());
    const size_t n_ep_block = n_ep_per_block + (b_idx < n_ep_rest ? 1 : 0);

    auto block_soln = eigensystem_hermitian(*ita, *itb, n_ep_block, map);
    assert_dbg(block_soln.n_ep() == n_ep_block, krims::ExcInternalError());

    // Lambda to pad the vector with zeros at beginning and end
    auto pad_with_zeros = [&begin_index, &A](const evector_type& v) {
      evector_type padded(A.size());
      std::copy(v.begin(), v.end(), padded.begin() + begin_index);
      return padded;
    };

    // Copy into solution data structure and pad the eigenvectors with zeros
    std::copy(block_soln.evalues().begin(), block_soln.evalues().end(),
              std::back_inserter(ret.evalues()));
    std::transform(block_soln.evectors().begin(), block_soln.evectors().end(),
                   std::back_inserter(ret.evectors()), pad_with_zeros);

    begin_index += ita->size();
  }
  assert_dbg(begin_index == A.size(), krims::ExcInternalError());

  return ret;
}

template <typename Matrix>
EigensolutionTypeFor<false, Matrix> eigensystem(
      const Matrix& A,
      typename std::enable_if<IsMatrix<Matrix>::value, typename Matrix::size_type>::type
            n_ep,
      const krims::GenMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, Matrix> problem_type;
  problem_type problem{A, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

template <typename MatrixA, typename MatrixB>
EigensolutionTypeFor<false, MatrixA> eigensystem(
      const MatrixA& A, const MatrixB& B,
      typename std::enable_if<IsMatrix<MatrixA>::value && IsMatrix<MatrixB>::value,
                              typename MatrixA::size_type>::type n_ep,
      const krims::GenMap& map) {
  // TODO This code is untested!
  assert_sufficiently_tested(false);

  // Setup problem
  typedef Eigenproblem<false, MatrixA, MatrixB> problem_type;
  problem_type problem{A, B, n_ep};
  return EigensystemSolver<problem_type>{map}.solve(problem).eigensolution();
}

}  // namespace linalgwrap
