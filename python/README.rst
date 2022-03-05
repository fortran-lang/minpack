Minpack Python bindings
=======================

Python bindings for Minpack.


Building the extension module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform a build some version of ``minpack`` has to be available on your system and preferably findable by ``pkg-config``.
Try to find a ``minpack`` installation you build against first with

.. code:: sh

   pkg-config --modversion minpack

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.


Using pip
^^^^^^^^^

This project support installation with pip as an easy way to build the Python API.

- C compiler to build the C-API and compile the extension module (the compiler name should be exported in the ``CC`` environment variable)
- Python 3.6 or newer
- The following Python packages are required additionally

  - `cffi <https://cffi.readthedocs.io/>`_
  - `numpy <https://numpy.org/>`_
  - `pkgconfig <https://pypi.org/project/pkgconfig/>`_ (setup only)

Make sure to have your C compiler set to the ``CC`` environment variable

.. code:: sh

   export CC=gcc

Install the project with pip

.. code:: sh

   pip install .

If you already have a ``minpack`` installation, *e.g.* from conda-forge, you can build the Python extension module directly without cloning this repository

.. code:: sh

   pip install "https://github.com/fortran-lang/minpack/archive/refs/heads/main.zip#egg=minpack&subdirectory=python"



Using meson
^^^^^^^^^^^

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
The out-of-tree build requires

- C compiler to build the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.53 or newer
- a build-system backend, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI <https://cffi.readthedocs.io/>`_ package installed

Setup a build with

.. code:: sh

   meson setup _build -Dpython_version=$(which python3)

The Python version can be used to select a different Python version, it defaults to ``'python3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build/minpack/_libminpack.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build
