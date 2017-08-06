//
// Copyright (C) 2016-17 by the lazyten authors
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

#pragma once
#include <functional>
#include <krims/NumComp.hh>
#include <rapidcheck.h>

namespace lazyten {
namespace tests {
using namespace krims;

/** \brief Helper class to simplify setup and calling comparative test
 * functions from rapidcheck.
 *
 * The idea is that we use rapidcheck to randomly create some test state,
 * represented by the type Args. Then we use generator functions
 * model_generator and sut_generator to generate an equivalent model
 * and system under test object.
 * These two are passed to comparative functions of the signature
 * ``testfunction_type``, which check whether the model and the sut
 * really *are* equivalent.
 *
 * The ``generate()`` function fuses the generator functions together
 * and returns a lambda, which lazily calls the test function.
 *
 * The ``run_test()`` function uses the lambda returned by ``generate()``
 * to actually perform a test.
 *
 * \tparam Args The argument type(s) we need to generate a sut matrix
 * or a model matrix (e.g. a matrix and a factor or a matrix
 * and some sizes ...). Use a std::tuple if more than one arg is needed for the
 * test state.
 */
template <typename Model, typename Sut, typename Args>
struct RCTestableGenerator {
  /** Type of the model implementation we test against */
  typedef Model model_type;

  /** Type of the system under test we want to test */
  typedef Sut sut_type;

  /** Type of the arguments, which are needed to construct a model and a sut
   * object */
  typedef Args args_type;

  /** Type of a testfunction template, which implements a comparative test. */
  typedef std::function<void(const model_type&, const sut_type&,
                             const NumCompAccuracyLevel)>
        testfunction_type;

  /** Type of a predicate function, which we can use to discard certain test
   *  cases if the args to generate the sut/model objects are not suitable
   *  for this test */
  typedef std::function<bool(const args_type&)> predfunction_type;

  /** \brief Construct a RapidcheckTestableGenerator object.
   *
   * \param sut_generator  A function that derives a System-under-test
   *                      matrix from an randomly generated core model
   *                      matrix of type compmat_type.
   * \param model_generator  A function that derives a compmat_type model
   *                        matrix in the same state as the test matrix
   *                        yielded by sutgenerator when applied to the
   *                        same core model matrix.
   */
  RCTestableGenerator(std::function<args_type(void)> args_generator_,
                      std::function<model_type(args_type)> model_generator_,
                      std::function<sut_type(args_type)> sut_generator_)
        : args_generator(args_generator_),
          model_generator(model_generator_),
          sut_generator(sut_generator_) {}

  /** \brief Construct a RCTestableGenerator object.
   *
   * Use explicit conversion to the sut_type and the model_type to obtain
   * the sut and the model from the args
   */
  RCTestableGenerator(std::function<args_type(void)> args_generator_)
        : args_generator(args_generator_),
          model_generator([](args_type t) { return model_type{t}; }),
          sut_generator([](args_type t) { return sut_type{t}; }) {}

  /** \brief Construct a RCTestableGenerator object.
   *
   * Use explicit conversion to the model_type to obtain
   * the model from the args
   */
  RCTestableGenerator(std::function<args_type(void)> args_generator_,
                      std::function<sut_type(args_type)> sut_generator_)
        : args_generator(args_generator_),
          model_generator([](args_type t) { return model_type{t}; }),
          sut_generator(sut_generator_) {}

  /* Return a std::function object that creates random arguments of type
   * args_type, then passes them to the model_generator and the sut_generator
   * to generate a model and a sut and then calles the testfunction func
   * in order to test whether the model and sut agree. r
   *
   * \param tolerance  The numeric tolerance level to use.
   * */
  std::function<void(void)> generate(
        testfunction_type func,
        const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default,
        predfunction_type pred = [](args_type) { return true; }) const;

  /** Generate a testable using ``generate()`` and run it using rapidcheck
   *
   * \param description for Rapidcheck
   * \param tolerance  The numeric tolerance level to use.
   * \param pred       Predicate which has to hold for the generated
   *                   arguments for the test case to be run.
   * */
  bool run_test(std::string description, testfunction_type func,
                const NumCompAccuracyLevel tolerance = NumCompAccuracyLevel::Default,
                predfunction_type pred = [](args_type) { return true; }) const {
    return rc::check(description, generate(func, tolerance, std::move(pred)));
  }

  std::function<args_type(void)> args_generator;
  std::function<model_type(args_type)> model_generator;
  std::function<sut_type(args_type)> sut_generator;
};

//
//
//
template <typename Model, typename Sut, typename Args>
std::function<void(void)> RCTestableGenerator<Model, Sut, Args>::generate(
      testfunction_type func, const NumCompAccuracyLevel tolerance,
      predfunction_type pred) const {
  // Return a lambda which calls the test appropriately.
  return [=] {
    // Generate the args
    Args arg = args_generator();

    // Run this case only if the predicate is true
    RC_PRE(pred(arg));

    // Call test function with appropriate model and sut:
    func(model_generator(arg), sut_generator(arg), tolerance);
  };
}

}  // namespace tests
}  // namespace lazyten
