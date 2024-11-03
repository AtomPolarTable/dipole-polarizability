.. _usr_installation:

Installation
############

Downloading the Code
====================

The latest code can be obtained from GitHub:

.. code-block:: bash

   git clone https://github.com/AtomPolarTable/dipole-polarizability.git

Installing
==========

To install `Pydipole`:

.. code-block:: bash

   pip install pydipole

To install the latest version from the repository:

.. code-block:: bash

   git clone https://github.com/AtomPolarTable/dipole-polarizability.git
   cd dipole-polarizability
   pip install .

Running Tests
=============

If you want to run the tests, install the test dependencies:

.. code-block:: bash

   pip install .[tests]

For developers who require all dependencies, install them with:

.. code-block:: bash

   pip install -e .[dev,tests]

Building the Documentation
==========================

To build the documentation using Sphinx into the `_build` directory:

.. code-block:: bash

    cd ./doc
    ./gen_api.sh
    sphinx-build -b html . _build

Dependencies
============

The following dependencies are required for `Pydipole` to build properly:

- `pandas` : https://github.com/pandas-dev/pandas
