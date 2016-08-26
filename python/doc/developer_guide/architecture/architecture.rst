Architecture considerations
===========================

Dependencies
------------

Several dependencies are needed in order to build the module:

 - OpenTURNS
 - Sphinx-doc (optional for this doc)

Compilation
-----------

.. code-block:: bash

    cd otsubsetinverse
    mkdir -p build && cd build
    cmake \
      -DCMAKE_INSTALL_PREFIX=$HOME/.local \
      -DOpenTURNS_DIR=<pythonpathlib>/lib/cmake/openturns \
      ..
    make
    make install
    make sphinx_pdf (optional for the pdf)
    make check
    make installcheck

