# linalgwrap
[![Build Status](https://travis-ci.org/linalgwrap/linalgwrap.svg?branch=master)](https://travis-ci.org/linalgwrap/linalgwrap) [![Licence](https://img.shields.io/github/license/linalgwrap/linalgwrap.svg)](LICENCE)

A lightweight linear algebra wrapper library adding support for lazy matrices
and lazy matrix evaluation to existing linear algebra libraries.

Note that this library is at a *very early stage* at the moment.
This means that interfaces can still change from the present situation
and only a fraction of the planned features is currently implemented.

## Dependencies
``linalgwrap`` depends on the following libraries:
- [krims](https://linalgwrap.org/krims) for many basic utilities
  (ParameterMap, Exception handling, Subscription pointers)
- [armadillo](http://arma.sourceforge.net/) as the only supported
  linear-algebra backend (so far)

Testing linalgwrap further requires
- [Catch](https://github.com/philsquared/Catch/) for the testing environment
- [rapidcheck](https://github.com/emil-e/rapidcheck) for property-based testing

Note, that for building ``linalgwrap`` (see below) you really only need to have
[armadillo](http://arma.sourceforge.net/) installed on your system.
All other dependencies can be automatically downloaded during the build process
if you choose to so (set ``AUTOCHECKOUT_MISSING_REPOS`` to ``ON``,
more below)

## Building
All compilers starting from ``clang-3.5`` and ``gcc-4.8`` should be able to build the code.
``C++11`` support is required and enables all basic functionality of the library.
Some things (most notably Block diagonal matrices with more than 4 blocks)
require ``C++14``, however.

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
ctest
```

## Short description of ``linalgwrap``
This section gives a very short description of the features of
linalgwrap some more detailed documentation will be produced.

Some of the design concepts and ideas that lead to the development
of linalgwrap can also be found in the material on
[blog.mfhs.eu](http://blog.mfhs.eu/uploads-publications/#Linalgwrap),
most notably the [MWM 2016 Poster](http://docs.mfhs.eu/conferences/2016_mwm/linalgwrap_lazy_linear_algebra_library.pdf).

### BaseInterfaces
- The classes in [src/linalgwrap/BaseInterfaces](BaseInterfaces)
  define the vector interface which is used inside linalgwrap.
- Fallback implementations for many important operations are provided
  in case a particular LA backend does not support these.
  This also minimises the effort to get a first working link to a new
  linear algebra backend.
- Implementations of backends (e.g. [src/linalgwrap/Armadillo](Armadillo))
  use this interface and specialise the fallback implementations
  such that the linear algebra backend performs the actual work.

### Builtin
- [src/linalgwrap/Builtin](Builtin) contains a builtin vector class,
  which is used as fallback.

### Lazy matrices
- Lazy matrices are a generalisation of "normal" matrices.
- The offer a matrix-like interface
  (addition, multiplication, application to vectors),
  but their elements may not neccessarily originate from
  a place in main memory.
- In other words lazy matrix elements may be computed by
  arbitrary computation on-the-fly while performing a
  matrix operation.
- The matrix may have state, which may be altered (updated).
- This (generally) makes obtaining individual matrix elements
  computationally less favourable than the application to
  a vector.
- Lazy matrices are subject to lazy (delayed) evaluation
  Lazy matrix expressions are only evaluated if a
  vector-apply is performed or if the
  user explicilty asks for it.
- Currently the [src/linalgwrap/DiagonalMatrix.hh](DiagonalMatrix)
  is available as an example of a builtin lazy matrix.
- Another example is available in the folder [examples/diagonal](examples/diagonal),
  namely the class [examples/diagonal/DiagonalUpdatable.hh](DiagonalUpdatable).
  This class also shows that custom lazy matrices can be created
  by simply inheriting from [src/linalgwrap/LazyMatrix_i.hh](LazyMatrix_i).

### TestingUtils
This class contains utilities for performing numerics-aware
property-based testing. This includes:
- An extension of ``krims::NumComp`` for comparing Matrices
  within error bounds (in file [src/linalgwrap/TestingUtils/krims_NumComp.hh](src/linalgwrap/TestingUtils/krims_NumComp.hh))
- Generators for values, vectors and matrices which are
  not too numerically difficult to deal with,
  such that tests do not fail due to numerics
  in folder [src/linalgwrap/TestingUtils/gen](src/linalgwrap/TestingUtils/gen))
