Install
============

bfieldtools is a pure Python package, and depends on several other Python packages. bfieldtools requires Python version 3.7 or higher.



Using pip
-----------

bfieldtools can most easily be installed using the pip package manager.

.. code-block:: bash

    pip install bfieldtools


A note on Python environments and dependencies
----------------------------------------------

Three-dimensional plots use PyVista, which is built on VTK. PyVista is installed
with bfieldtools. For interactive plots in Jupyter, install the notebook extras:

.. code-block:: bash

    pip install "pyvista[jupyter]>=0.48,<0.49"

Select the trame backend in a notebook:

.. code-block:: python

    import pyvista as pv
    pv.set_jupyter_backend("trame")

The documentation includes static images and interactive scenes. The interactive
scenes run in the browser without a Python server; Python callbacks and fresh
calculations require a running notebook.

To build the documentation locally from this checkout, use Python 3.10 or
newer (the uv configuration defaults to Python 3.11):

.. code-block:: bash

    uv sync --locked --extra docs
    uv run --locked --extra docs make -C docs html

The output is in ``docs/_build/html``. The documentation build renders plots
off screen and needs a working VTK/OpenGL rendering environment. External-data
examples, the numerical benchmarks, and large solver-dependent optimizations are
included as source-only pages in the default build. Set
``BFIELDTOOLS_DOCS_FULL_GALLERY=1`` to run the large optimizations using
open-source solvers; see ``source_only_examples`` in ``docs/conf.py``.

From source
-----------

To install the *development* version of bfieldtools, run

.. code-block:: bash

    wget https://github.com/bfieldtools/bfieldtools/archive/master.zip

After downloading and unzipping the archive, basic setuptools workflow applies.

.. code-block:: bash

    cd bfieldtools/
    python setup.py install


Additional solvers
-------------------

The stream function optimization functionality of bfieldtools relies on external numerical solvers, for which the installation procedures may vary.  For more details, see the `documentation of CVXPY`_. OSQP, CVXOPT, SCS and ECOS are installed with bfieldtools and are also available through conda and pip. MOSEK_ is (optional solver dependency, recommended) robust and fast, but commercial (requires license). Installable via pip or conda.

.. _documentation of CVXPY: https://www.cvxpy.org/install/index.html#install-from-source

.. _MOSEK: https://docs.mosek.com/9.0/pythonapi/install-interface.html
