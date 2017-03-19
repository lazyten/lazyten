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

#pragma once
#include <complex>
#include <cstddef>
#include <iterator>
#include <krims/Algorithm.hh>
#include <krims/ExceptionSystem.hh>
#include <krims/TypeUtils.hh>

namespace linalgwrap {
DefException1(ExcUnknownEvalWhich, std::string, << "The vaule " << arg1
                                                << " for which is not known");

/** Comparator to use to canonically sort eigenvalues */
struct EvalComp {
  template <typename T, bool iscpx = krims::IsComplexNumber<T>::value>
  bool operator()(typename std::enable_if<iscpx, T>::type i, T j) const {
    return std::abs(i) < std::abs(j);
  }

  template <typename T, bool iscpx = krims::IsComplexNumber<T>::value>
  bool operator()(typename std::enable_if<!iscpx, T>::type i, T j) const {
    return i < j;
  }
};

/** Select the n_ep eigenvalues matching the which criterion and
 * return the index array which would produce them in sorted manor
 * if drawn from the container.
 *
 * The sorting for real eigenvalues is by value and for complex ones
 * by magnitude.
 *
 * This function is mainly intended for dealing with direct eigensolver
 * methods which do not always return a sorted list of eigenpairs and
 * moreover sometimes return many more eigenvalues than desired.
 *
 * \param eval_begin   The beginning of the eigenvalue range to consider
 * \param eval_end     The end of the eigenvalue range to consider
 * \param which        which eigenpairs to keep
 * \param n_ep         Number of eigenpairs to keep
 * \param sorted_canonically   May the algorithm assume that the eigenvalues
 *                     are already sorted canonically, i.e. in the way
 *                     the EvalComp functor above would sort them.
 */
template <typename Iterator>
std::vector<size_t> select_eigenvalues(Iterator eval_begin, Iterator eval_end,
                                       const std::string& which, const size_t n_ep,
                                       const bool sorted_canonically = false) {
  assert_dbg(which == "SM" || which == "LM" || which == "LR" || which == "SR" ||
                   which == "SI" || which == "LI",
             ExcUnknownEvalWhich(which));

  typedef typename std::iterator_traits<Iterator>::value_type evalue_type;

  // The size of the eigenvalue range
  const size_t size = static_cast<size_t>(std::distance(eval_begin, eval_end));

  // Is the eigenvalue type a real type
  constexpr bool real = !krims::IsComplexNumber<evalue_type>::value;

  // Is canonical sorting requested:
  const bool want_canonical = (real && which[1] == 'R') || (!real && which[1] == 'M');

  // Indices of the eigenpairs which should be kept,
  // given in the order they should appear
  std::vector<size_t> idcs;

  if (want_canonical && sorted_canonically) {
    // We want canonical sorting and the evals are already sorted
    // => no work to do
    idcs.resize(size);
    std::iota(std::begin(idcs), std::end(idcs), 0);
  } else {
    using std::abs;
    using std::real;
    using std::imag;

    if (which[1] == 'R') {
      // Sort by real part
      auto comp = [](const evalue_type& lhs, const evalue_type& rhs) {
        return real(lhs) < real(rhs);
      };
      idcs = krims::argsort(eval_begin, eval_end, std::move(comp));
    } else if (which[1] == 'I') {
      // Sort by imag. part
      auto comp = [](const evalue_type& lhs, const evalue_type& rhs) {
        return imag(lhs) < imag(rhs);
      };
      idcs = krims::argsort(eval_begin, eval_end, std::move(comp));
    } else {
      // Sort by magnitude
      auto comp = [](const evalue_type& lhs, const evalue_type& rhs) {
        return abs(lhs) < abs(rhs);
      };
      idcs = krims::argsort(eval_begin, eval_end, std::move(comp));
    }
  }

  // Select the eigenvalues we care about:
  // If we want to keep the smallest we need the beginning
  // of the range, else the end
  size_t start = which[0] == 'S' ? 0 : size - n_ep;
  size_t end = which[0] == 'S' ? n_ep : size;

  if (want_canonical) {
    // The ordering we just did above is the canonical ordering,
    // which EvalComp would do anyway:
    // => just shrink the vector to the entries we care about
    idcs = std::vector<size_t>(std::begin(idcs) + static_cast<ptrdiff_t>(start),
                               std::begin(idcs) + static_cast<ptrdiff_t>(end));
  } else {
    // Build a temporary evalues array:
    std::vector<evalue_type> tmp;
    tmp.reserve(end - start);
    for (size_t i = start; i < end; ++i) {
      tmp.push_back(*(eval_begin + idcs[i]));
    }

    // Argsort it using the default eigenpair ordering
    std::vector<size_t> i2 = krims::argsort(std::begin(tmp), std::end(tmp), EvalComp{});

    // Translate such that the desired order is included
    std::vector<size_t> idcs_new;
    idcs_new.reserve(i2.size());
    for (const auto& i : i2) {
      idcs_new.push_back(idcs[start + i]);
    }
    idcs = idcs_new;
  }
  assert_dbg(idcs.size() == n_ep, krims::ExcInternalError());
  return idcs;
}

}  // namespace linalgwrap
