# lazyten
[![Build Status](https://travis-ci.org/lazyten/lazyten.svg?branch=master)](https://travis-ci.org/lazyten/lazyten)
[![Coverage Status](https://coveralls.io/repos/github/lazyten/lazyten/badge.svg?branch=master)](https://coveralls.io/github/lazyten/lazyten)
[![Licence](https://img.shields.io/github/license/lazyten/lazyten.svg)](LICENCE)

A lightweight linear algebra wrapper library adding support for
[lazy matrices](#lazy-matrices)
and lazy matrix evaluation to existing linear algebra libraries.
This library used to be called `linalgwrap` and was just recently renamed to `lazyten`
(for **lazy** **ten**sor library).
Right now we are only able to deal with lazy matrices,
but we intend to add support for lazy tensor evaluation as well in the future.

Note that we are still at a *very early* development stage.
In other words interfaces will very likely change in the future
in incompatible ways to what is currently implemented.
We try to make sure that this does, however,
not go unnoticed, i.e. that existing code breaks loudly at compile time.

## Dependencies
``lazyten`` depends on the following libraries:
- [krims](https://lazyten.org/krims) for many basic utilities
  (GenMap, Exception handling, Subscription pointers)
- A BLAS implementation, e.g. [OpenBLAS](https://github.com/xianyi/OpenBLAS/)
- A LAPACK compatible library, e.g.
  [LAPACK](http://netlib.org/lapack) in order to use the LAPACK eigensolvers
- [armadillo](http://arma.sourceforge.net/) for the armadillo eigensolvers
  and linear solvers as well as the only linear-algebra backend (so far)
- *(optional)* [ARPACK](http://www.caam.rice.edu/software/ARPACK/) in order to use
  ``lazyten`` with the ARPACK eigensolver.

Testing ``lazyten`` further requires
- [Catch](https://github.com/philsquared/Catch/) for the testing environment
- [rapidcheck](https://github.com/emil-e/rapidcheck) for property-based testing

`lazyten` comes with a couple of tools to automatically download and build
some of these dependencies during the
[manual build process](#building-manually-without-spack).
To avoid the potential hassle all together we do, however,
recommend [using `Spack` to build `lazyten`](#building-via-spack-recommended).


## Building via Spack (recommended)
The [Spack](https://spack.io) package manager allows to easily install scientific software.
Both `lazyten` as well as `krims` are available in `spack`.
Installing `lazyten` via `spack` is therefore as simple as
```sh
# Clone and setup spack
git clone https://github.com/llnl/spack.git
export SPACK_ROOT="$PWD/spack"
. $SPACK_ROOT/share/spack/setup-env.sh

# Install lazyten (and all required dependencies)
spack install lazyten
```
Once this has happened you can add all relevant environment variables
(`LD_LIBRARY_PATH`, `PATH`, `CPATH`, ...) to the current shell
via the commands
```sh
spack module loads -r lazyten > /tmp/lazyten.modules
. /tmp/lazyten.modules
```
which will generate a list of all spack modules `lazyten` needs
and loades them thereafter.
Running the above two lines of code gets you ready for
linking `lazyten` to your project.

Other than that Spack makes it very easy to customise the installation, too.
For example to influence which features of `lazyten` are to be built,
one can add specifiers to the `spack install` command.  
If you prefer to build `lazyten` without `ARPACK` for example, run
```sh
spack install lazyten~arpack
```
instead of the command mentioned initially.

For more information about how to use Spack see
the great [Spack documentation](https://spack.readthedocs.io).


## Building manually (without Spack)
For configuring the build we need at least ``cmake`` ``3.0.0``.  
All compilers starting from ``clang-3.5`` and ``gcc-4.8`` should be able to build the code.
``C++11`` support is required and enables all basic functionality of the library.
A couple of things require ``C++14``, however.

If you choose to build with the flag ``AUTOCHECKOUT_MISSING_REPOS`` set to ``ON``
most required dependencies (*except* armadillo, LAPACK and BLAS) will be automatically
downloaded and compiled alongside ``lazyten``,
so you really need to have only armadillo, LAPACK and BLAS on your system at build time.

In order to build ``lazyten`` with tests (recommended) run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON ..
cmake --build .
ctest
```

In order to build without tests run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON -DLAZYTEN_ENABLE_TESTS=OFF -DKRIMS_ENABLE_TESTS=OFF ..
cmake --build .
```

## Short description of ``lazyten``
This section gives a very short description of the features of
``lazyten``.
We hope to produce some more detailed documentation at some point.

Some of the design concepts and ideas that lead to the development
of the present version of ``lazyten`` can also be found in the material on
[michael-herbst.com](https://michael-herbst.com/tag/lazy-matrices.html),
most notably the presentation at the Niels Bohr Institute
[HPC Day 2017](https://michael-herbst.com/talks/2017.05.19_HPC_Day_NBI.pdf)
as well as parts of the [Design of ``molsturm``](http://docs.mfhs.eu/phd/invited_talks/2016.12.09_Design_Molsturm.pdf)
talk.

### BaseInterfaces
- The classes in [BaseInterfaces](src/lazyten/BaseInterfaces)
  define the vector interface which is used inside ``lazyten``.
- Fallback implementations for many important operations are provided
  in case a particular linear-algebra backend does not support these.
  This also minimises the effort to get a first working link to a new
  linear algebra backend.
- Implementations of backends (e.g. [Armadillo](src/lazyten/Armadillo))
  use this interface and specialise the default implementations
  such that the backend library performs the actual work.
- Our goal is to forward as much of the performance optimisations of the
  backend libraries to the interface, without depending on only
  one backend.

### Builtin
- [Builtin](src/lazyten/Builtin) contains a builtin vector class,
  which is used as fallback.

### Lazy matrices
- Lazy matrices are a generalisation of "normal" matrices.
- They offer a matrix-like interface
  (addition, multiplication, application to vectors),
  but their elements do not need to be placed in a
  stride of memory.
- In other words lazy matrix elements may be computed by
  arbitrary computation on-the-fly while performing a
  matrix operation.
- The matrix may have state, which may be altered (updated).
- This (generally) makes obtaining individual matrix elements
  computationally less favourable than the application to
  a vector.
- Lazy matrices are subject to lazy (delayed) evaluation.
  Lazy matrix expressions are only evaluated if a
  vector-apply is performed or if the
  user explicitly asks for it.
- Currently the [DiagonalMatrix](src/lazyten/DiagonalMatrix.hh)
  is available as an example of a builtin lazy matrix.
- Another example is the class [DiagonalUpdatable](examples/diagonal/DiagonalUpdatable.hh)
  of the [diagonal](examples/diagonal) example program.
  This class also shows that custom lazy matrices can be created
  by simply inheriting from [LazyMatrix_i](src/lazyten/LazyMatrix_i.hh).

### Matrix operations
#### As Matrix member functions
All matrices support the following member functions:
- ``apply``: Generalised matrix-vector product (like the ``gemv`` BLAS call).
  Allows to multiply a matrix with a number of vectors and add (or set) the result
  to a different set of vectors.
- ``mmult``: Generalised matrix-matrix product (similar to BLAS' ``gemm``).
  Allows to multiply two matrices and add or set the result to a third.
- ``extract_block``: Extract a block of matrix values and add or set them
  to some pre-allocated storage

#### Matrix operators and global scope
On global scope we have:
- ``as_stored(mat)``: Either return a reference to the current object (in case it already is a stored matrix)
  or return a stored representation (i.e. a copy) of the matrix ``mat``.
- ``trans(mat)``: Return an object which represents the transpose of a matrix.
- ``conjtrans(mat)``: Return an object which represents the conjugate transpose of a matrix.
- ``operator*``: Multiplication between matrices and matrix and vector.
  Internally calls the ``apply`` and ``mmult`` methods above.
- All kinds of matrix norms: ``norm_linf``, ``norm_l1``, ``norm_frobenius``, ``norm_frobenius_squared``.

### Solvers
- In order to solve an eigenproblem the methods ``eigensystem`` and ``eigensystem_hermitian``
  in the file [eigensystem.hh](src/lazyten/eigensystem.hh) are available as
  high-level routines. These are the recommended routines to solve eigenproblems,
  since their interface is designed to be easy to use and it easy to enforce
  a particular eigensolver explicitly (using the parameter key ``method``).
- For linear problems the file [solve.hh](src/lazyten/solve.hh) similarly
  contains the method  ``solve``.
- The folder [Base/Solvers](src/lazyten/Base/Solvers) holds all the lower
  interfaces and some utilities the actual solvers use.
- Currently only [ArpackEigensolver](src/lazyten/Arpack/ArpackEigensolver.hh)
  and some eigensolvers from Lapack (either indirectly via
  [ArmadilloEigensolver](src/lazyten/Armadillo/ArmadilloEigensolver.hh)
  or directy via [LapackEigensolver](src/lazyten/Lapack/LapackEigensolver.hh))
  are implemented as eigensolver backends.
- ARPACK is enabled if it is found on the system.
- Linear problems are always solved with ``LAPACK`` via  ``armadillo`` right now.

### TestingUtils
This class contains utilities for performing numerics-aware
property-based testing. This includes:
- An extension of [``krims::NumComp``](https://lazyten.org/krims/#performing-floating-point-comparisons)
  for comparing Matrices within error bounds (in file [TestingUtils/krims_NumComp.hh](src/lazyten/TestingUtils/krims_NumComp.hh))
- Generators for scalar values, vectors and matrices which are
  not too difficult to deal with,
  such that tests do not fail due to accumulation of numeric errors
  (in folder [TestingUtils/gen](src/lazyten/TestingUtils/gen))
