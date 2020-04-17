Install
============

bfieldtools is a pure Python package, and depends on several other Python packages. bfieldtools requires Python version 3.7 or higher.  


Using conda
-----------

We strongly recommend the Anaconda_ distribution of Python, which comes with more than 250 scientific packages pre-bundled, and includes the conda command line tool for installing new packages and managing different package sets (“environments”) for different projects. Alternatively, the more bare-bones Miniconda_ distribution should also work.


.. _Anaconda: https://www.anaconda.com/
.. _Miniconda: : https://docs.conda.io/en/latest/miniconda.html

.. code-block:: bash

    WIP


Using pip
-----------

bfieldtools can also be installed using the pip package manager. 

.. code-block:: bash

    pip install bfieldtools

From source
-----------

Basic setuptools workflow applies.
    
.. code-block:: bash

    python setup.py install


Additional solvers
-------------------

The stream function optimization functionality of bfieldtools relies on external numerical solvers, for which the installation procedures may vary. Here, we link to the installation procedures of these solvers.

MOSEK_
^^^^^^^^
Optional solver dependecy. Robust and fast, but commercial (requires licenese). Installable via pip or conda.

OSQP
^^^^^^
Installed with bfieldtools, available through conda and pip

CVXOPT
^^^^^^
Installed with bfieldtools, available through conda and pip

SCS
^^^^^^
Installed with bfieldtools, available through conda and pip
 
ECOS
^^^^^^
Installed with bfieldtools, available through conda and pip


.. _MOSEK: https://docs.mosek.com/9.0/pythonapi/install-interface.html
