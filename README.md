# linalgwrap
[![Build Status](https://travis-ci.org/linalgwrap/linalgwrap.svg?branch=master)](https://travis-ci.org/linalgwrap/linalgwrap)
[![Coverage Status](https://coveralls.io/repos/github/linalgwrap/linalgwrap/badge.svg?branch=master)](https://coveralls.io/github/linalgwrap/linalgwrap)
[![Licence](https://img.shields.io/github/license/linalgwrap/linalgwrap.svg)](LICENCE)

A lightweight linear algebra wrapper library adding support for [lazy-matrices](#lazy matrices)
and lazy matrix evaluation to existing linear algebra libraries.

Note that this library is at a *very early stage* at the moment.
This means that interfaces will very likely change in the future
and only a fraction of the planned features are currently implemented.

## Dependencies
``linalgwrap`` depends on the following libraries:
- [krims](https://linalgwrap.org/krims) for many basic utilities
  (GenMap, Exception handling, Subscription pointers)
- A BLAS implementation, e.g. [OpenBLAS](https://github.com/xianyi/OpenBLAS/)
- A LAPACK compatible library, e.g.
  [LAPACK](http://netlib.org/lapack) in order to use the LAPACK eigensolvers
- [armadillo](http://arma.sourceforge.net/) for the armadillo eigensolvers
  and linear solvers as well as the only linear-algebra backend (so far)
- *(optional)* [ARPACK](http://www.caam.rice.edu/software/ARPACK/) in order to use
  ``linalgwrap`` with the ARPACK eigensolver.

Testing ``linalgwrap`` further requires
- [Catch](https://github.com/philsquared/Catch/) for the testing environment
- [rapidcheck](https://github.com/emil-e/rapidcheck) for property-based testing

Note, that for building ``linalgwrap`` (see [below](#building)) you really only need to have
[armadillo](http://arma.sourceforge.net/), [LAPACK](http://netlib.org/lapack) and a BLAS
installed on your system.
All other dependencies can be automatically downloaded during the build process
if you choose to do so (set ``AUTOCHECKOUT_MISSING_REPOS`` to ``ON``,
more below)

## Building
For configuring the build we need at least ``cmake`` ``3.0.0``.  
All compilers starting from ``clang-3.5`` and ``gcc-4.8`` should be able to build the code.
``C++11`` support is required and enables all basic functionality of the library.
A couple of things require ``C++14``, however.

If you choose to build with the flag ``AUTOCHECKOUT_MISSING_REPOS`` set to ``ON``
all required dependencies (**except** armadillo) will be automatically downloaded
and compiled alongside ``linalgwrap``.

In order to build ``linalgwrap`` with tests (recommended) run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON ..
cmake --build .
ctest
```

In order to build without tests run
```
mkdir build && cd build
cmake -DAUTOCHECKOUT_MISSING_REPOS=ON -DLINALGWRAP_ENABLE_TESTS=OFF -DKRIMS_ENABLE_TESTS=OFF ..
cmake --build .
```

## Short description of ``linalgwrap``
This section gives a very short description of the features of
``linalgwrap``.
We hope to produce some more detailed documentation at some point.

Some of the design concepts and ideas that lead to the development
of the present version of ``linalgwrap`` can also be found in the material on
[michael-herbst.com](https://michael-herbst.com/tag/lazy-matrices.html),
most notably the presentation at the Niels Bohr Institute
[HPC Day 2017](https://michael-herbst.com/talks/2017.05.19_HPC_Day_NBI.pdf)
as well as parts of the [Design of ``molsturm``](http://docs.mfhs.eu/phd/invited_talks/2016.12.09_Design_Molsturm.pdf)
talk.

### BaseInterfaces
- The classes in [BaseInterfaces](src/linalgwrap/BaseInterfaces)
  define the vector interface which is used inside ``linalgwrap``.
- Fallback implementations for many important operations are provided
  in case a particular linear-algebra backend does not support these.
  This also minimises the effort to get a first working link to a new
  linear algebra backend.
- Implementations of backends (e.g. [Armadillo](src/linalgwrap/Armadillo))
  use this interface and specialise the default implementations
  such that the backend library performs the actual work.
- Our goal is to forward as much of the performance optimisations of the
  backend libraries to the interface, without depending on only
  one backend.

### Builtin
- [Builtin](src/linalgwrap/Builtin) contains a builtin vector class,
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
- Currently the [DiagonalMatrix](src/linalgwrap/DiagonalMatrix.hh)
  is available as an example of a builtin lazy matrix.
- Another example is the class [DiagonalUpdatable](examples/diagonal/DiagonalUpdatable.hh)
  of the [diagonal](examples/diagonal) example program.
  This class also shows that custom lazy matrices can be created
  by simply inheriting from [LazyMatrix_i](src/linalgwrap/LazyMatrix_i.hh).

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
  in the file [eigensystem.hh](src/linalgwrap/eigensystem.hh) are available as
  high-level routines. These are the recommended routines to solve eigenproblems,
  since their interface is designed to be easy to use and it easy to enforce
  a particular eigensolver explicitly (using the parameter key ``method``).
- For linear problems the file [solve.hh](src/linalgwrap/solve.hh) similarly
  contains the method  ``solve``.
- The folder [Base/Solvers](src/linalgwrap/Base/Solvers) holds all the lower
  interfaces and some utilities the actual solvers use.
- Currently only [ArpackEigensolver](src/linalgwrap/Arpack/ArpackEigensolver.hh)
  and some eigensolvers from Lapack (either indirectly via
  [ArmadilloEigensolver](src/linalgwrap/Armadillo/ArmadilloEigensolver.hh)
  or directy via [LapackEigensolver](src/linalgwrap/Lapack/LapackEigensolver.hh))
  are implemented as eigensolver backends.
- ARPACK is enabled if it is found on the system.
- Linear problems are always solved with ``LAPACK`` via  ``armadillo`` right now.

### TestingUtils
This class contains utilities for performing numerics-aware
property-based testing. This includes:
- An extension of [``krims::NumComp``](https://linalgwrap.org/krims/#performing-floating-point-comparisons)
  for comparing Matrices within error bounds (in file [TestingUtils/krims_NumComp.hh](src/linalgwrap/TestingUtils/krims_NumComp.hh))
- Generators for scalar values, vectors and matrices which are
  not too difficult to deal with,
  such that tests do not fail due to accumulation of numeric errors
  (in folder [TestingUtils/gen](src/linalgwrap/TestingUtils/gen))
