.. image:: https://travis-ci.org/openturns/otsubsetinverse.svg?branch=master
    :target: https://travis-ci.org/openturns/otsubsetinverse

OTSubsetInverse Module
======================

The goal of the module is to compute the performance of a function (threshold value) for a given target probability. The technique uses the subset sampling in the same way as in the original.

Compilation
===========

Linux
-----

    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=$PWD/install \
    -DOpenTURNS_DIR=OT_PREFIX/lib/cmake/openturns ..
    make
    make check
    make install
    make installcheck

Windows
-------

Use argument make PYBASEVER=3.6 for python 3.6.

- 32 bits

    cd otsubsetinverse/distro/windows
    make OT_PREFIX=..openturns-src/build-i686-w64-mingw32/install

- 64 bits

    make ARCH=x86_64 OT_PREFIX=..openturns-src/build-x86_64-w64-mingw32/install  # (for a 64 bits target)

