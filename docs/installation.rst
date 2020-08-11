Install
============

bfieldtools is a pure Python package, and depends on several other Python packages. bfieldtools requires Python version 3.6 or higher.



Using pip
-----------

bfieldtools can most easily be installed using the pip package manager. 

.. code-block:: bash

    pip install bfieldtools
    
    
A note on Python environments and dependencies
----------------------------------------------

bfieldtools requires some dependencies to work. The most demanding dependency is mayavi, which is used for 3D plots. Mayavi, which in turn relies on vtk, can be a hassle to install on some platforms. One of the best ways to getn bfieldtools fully working is to use the Anaconda_ or Miniconda_ Python distributions. You can then run 

.. code-block:: bash

    conda install mayavi
    
to install mayavi, after which you can install bfieldtools using pip as seen above. bfieldtools itself is unfortunately not installable through conda since some other dependencies are pip-only.


.. _Anaconda: https://www.anaconda.com/products/individual

.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html

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
