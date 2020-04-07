.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_validation_validate_line_charge_potential.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_validation_validate_line_charge_potential.py:


Line charge
================

Test and validation of gamma_0

For the math see:
        A. S. Ferguson, Xu Zhang and G. Stroink,
        "A complete linear discretization for calculating the magnetic field
        using the boundary element method,"
        in IEEE Transactions on Biomedical Engineering,
        vol. 41, no. 5, pp. 455-460, May 1994.
        doi: 10.1109/10.293220


.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    from mayavi import mlab

    path = "/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools"
    if path not in sys.path:
        sys.path.insert(0, path)

    from bfieldtools.integrals import gamma0
    from bfieldtools.utils import tri_normals_and_areas








%% Test potential shape slightly above the surface


.. code-block:: default

    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    tris = np.array([[0, 1, 2]])
    p_tris = points[tris]

    # Evaluation points
    Nx = 100
    xx = np.linspace(-2, 2, Nx)
    X, Y = np.meshgrid(xx, xx, indexing="ij")
    Z = np.zeros_like(X) + 0.1
    p_eval = np.array([X, Y, Z]).reshape(3, -1).T

    # Difference vectors
    RR = p_eval[:, None, None, :] - p_tris[None, :, :, :]

    pot = gamma0(RR)

    # Plot shape
    plt.figure()
    plt.imshow(
        pot[:, 0].sum(axis=-1).reshape(Nx, Nx),
        extent=(xx.min(), xx.max(), xx.max(), xx.min()),
    )
    plt.ylabel("x")
    plt.xlabel("y")




.. image:: /auto_examples/validation/images/sphx_glr_validate_line_charge_potential_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Text(0.5, 0, 'y')




.. code-block:: default

    """ Test potential at directly at the edge. As the line has no
        perpendicular dimensions, the potential is infinite at the edge.
        The regularization factor given to the function apprximates the line
        current with a small radius, giving rougly constant potential on the line.
        The relative error between the regularized verison and
        the infinitely thin line charge seems to be on the order of the "reg" value

        The "symmetrize" option symmeterizes the result with respect to the
        mid point. This removes errors (Nans) on the other continuoation of the edge
    """




.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    ' Test potential at directly at the edge. As the line has no\n    perpendicular dimensions, the potential is infinite at the edge.\n    The regularization factor given to the function apprximates the line\n    current with a small radius, giving rougly constant potential on the line.\n    The relative error between the regularized verison and\n    the infinitely thin line charge seems to be on the order of the "reg" value\n\n    The "symmetrize" option symmeterizes the result with respect to the\n    mid point. This removes errors (Nans) on the other continuoation of the edge\n'




.. code-block:: default

    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])

    tris = np.array([[0, 1, 2]])
    p_tris = points[tris]

    # Evaluation points
    Nx = 1000
    x = np.linspace(-2, 2, Nx)
    y = z = np.zeros_like(x)
    p_eval = np.array([x, y, z]).T

    # Difference vectors
    RR = p_eval[:, None, None, :] - p_tris[None, :, :, :]

    # Regularize and symmetrize
    pot = gamma0(RR, 1e-13, True)
    pot0 = pot[:, 0, 2]
    plt.figure()
    plt.plot(x, pot0, linewidth=5)
    # Symmetrize, but do not regularize
    pot = gamma0(RR, 0, True)
    pot1 = pot[:, 0, 2]
    plt.plot(x, pot1, "--", linewidth=3)
    # Neither
    pot = gamma0(RR, 0, False)
    pot2 = pot[:, 0, 2]
    plt.plot(x, pot2)
    plt.xlabel("x")
    plt.legend(("reg + sym", "sym", "neither"))

    plt.figure()
    plt.title("Relative error")
    plt.plot(x, abs(pot0 - pot1) / pot1)



.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/validation/images/sphx_glr_validate_line_charge_potential_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/validation/images/sphx_glr_validate_line_charge_potential_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\bfieldtools\integrals.py:73: RuntimeWarning: invalid value encountered in true_divide
      res = np.log((nn1 + dotprods1 + reg) / (nn2 + dotprods2 + reg))
    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\bfieldtools\integrals.py:73: RuntimeWarning: divide by zero encountered in log
      res = np.log((nn1 + dotprods1 + reg) / (nn2 + dotprods2 + reg))
    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\bfieldtools\integrals.py:81: RuntimeWarning: divide by zero encountered in true_divide
      (nn1[mask] - dotprods1[mask] + reg) / (nn2[mask] - dotprods2[mask] + reg)
    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\bfieldtools\integrals.py:73: RuntimeWarning: invalid value encountered in true_divide
      res = np.log((nn1 + dotprods1 + reg) / (nn2 + dotprods2 + reg))
    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\bfieldtools\integrals.py:73: RuntimeWarning: divide by zero encountered in log
      res = np.log((nn1 + dotprods1 + reg) / (nn2 + dotprods2 + reg))
    C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\examples\validation\validate_line_charge_potential.py:101: RuntimeWarning: invalid value encountered in true_divide
      plt.plot(x, abs(pot0 - pot1) / pot1)

    [<matplotlib.lines.Line2D object at 0x0000025409DA30B8>]




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.260 seconds)


.. _sphx_glr_download_auto_examples_validation_validate_line_charge_potential.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: validate_line_charge_potential.py <validate_line_charge_potential.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: validate_line_charge_potential.ipynb <validate_line_charge_potential.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
