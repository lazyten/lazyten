# linalgwrap
[![Build Status](https://travis-ci.org/linalgwrap/linalgwrap.svg?branch=master)](https://travis-ci.org/linalgwrap/linalgwrap) [![Licence](https://img.shields.io/github/license/linalgwrap/linalgwrap.svg)](LICENCE)

A lightweight linear algebra wrapper library adding support for lazy matrices
and lazy matrix evaluation to existing linear algebra libraries.

Note that this library is at a *very early stage* at the moment.
This means that interfaces can still change from the present situation.

## Dependencies
``linalgwrap`` depends on the following libraries:
- [krims](https://linalgwrap.org/krims) for many basic utilities
  (ParameterMap, Exception handling, Subscription pointers)
- [armadillo](http://arma.sourceforge.net/) as the linear-algebra backend

Testing linalgwrap further requires
- [Catch](https://github.com/philsquared/Catch/) for the testing environment
- [rapidcheck](https://github.com/emil-e/rapidcheck) for property-based testing

Note, that for building ``linalgwrap`` (see below) you really only need to have
[armadillo](http://arma.sourceforge.net/) installed on your system.
All other dependencies can be automatically downloaded during the build process
if you choose to so (set ``AUTOCHECKOUT_MISSING_REPOS`` to ``ON``,
more below)

## Building and dependencies
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
