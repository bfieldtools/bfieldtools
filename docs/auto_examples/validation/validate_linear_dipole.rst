.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_validation_validate_linear_dipole.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_validation_validate_linear_dipole.py:


Linear dipole
=============

Test and validation of potential of linearly distributed dipolar density

For the math see:
        J. C. de Munck, "A linear discretization of the volume conductor
        boundary integral equation using analytically integrated elements
        (electrophysiology application),"
        in IEEE Transactions on Biomedical Engineering,
        vol. 39, no. 9, pp. 986-990, Sept. 1992.
        doi: 10.1109/10.256433




.. code-block:: default


    import numpy as np
    import matplotlib.pyplot as plt
    import sys
    path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
    if path not in sys.path:
        sys.path.insert(0,path)

    from bfieldtools.integrals import triangle_potential_dipole_linear
    from bfieldtools.integrals import omega
    from bfieldtools.utils import tri_normals_and_areas







%% Test potential shape slightly above the surface
points = np.array([[0,0,0],
                   [1,0,0],
                   [0,1,0]])

tris = np.array([[0,1,2]])
p_tris = points[tris]


.. code-block:: default


    points = np.array([[0,0,0],
                       [1,1,0],
                       [1,-1,0],
                       [-1,-1,0],
                       [-1,1,0]])

    tris = np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,1]])
    tris = np.flip(tris,axis=-1)
    p_tris = points[tris]

    # Evaluation points
    Nx = 100
    xx = np.linspace(-2, 2, Nx)
    X,Y = np.meshgrid(xx, xx, indexing='ij')
    Z = np.zeros_like(X) + 0.01
    p_eval = np.array([X,Y,Z]).reshape(3,-1).T

    # Difference vectors
    RR = p_eval[:,None,None,:] - p_tris[None,:,:,:]
    tn, ta = tri_normals_and_areas(points, tris)

    pot = triangle_potential_dipole_linear(RR, tn, ta, False)

    # Plot shapes
    f, ax = plt.subplots(1, 3)
    for i in range(3):
        plt.sca(ax[i])
        plt.imshow(pot[:,2,i].reshape(Nx, Nx), extent=(xx.min(),xx.max(),
                                                       xx.max(),xx.min()))
        plt.colorbar(orientation='horizontal')
        if i==0:
            plt.ylabel('x')
            plt.xlabel('y')




.. image:: /auto_examples/validation/images/sphx_glr_validate_linear_dipole_001.png
    :class: sphx-glr-single-img




%% Test summation formula


.. code-block:: default

    pot_sum = triangle_potential_dipole_linear(RR, tn, ta, False).sum(axis=-1)
    solid_angle = omega(RR)

    # Plot shapes
    f, ax = plt.subplots(1, 3)
    plt.sca(ax[0])
    plt.title('Sum of potentials')
    plt.imshow(pot_sum[:,0].reshape(Nx, Nx), vmin=0, vmax=pot_sum.max())
    plt.colorbar(orientation='horizontal')
    plt.sca(ax[1])
    plt.title('Solid angle')
    plt.imshow(solid_angle[:,0].reshape(Nx, Nx), vmin=0, vmax=pot_sum.max())
    plt.colorbar(orientation='horizontal')
    plt.sca(ax[2])
    plt.title('Abs difference')
    plt.imshow(abs((solid_angle[:,0]-pot_sum[:,0])).reshape(Nx, Nx),
               vmin=0, vmax=pot_sum.max()/1e16)
    plt.colorbar(orientation='horizontal', pad=-0.2)
    plt.axis('image')

    plt.tight_layout()





.. image:: /auto_examples/validation/images/sphx_glr_validate_linear_dipole_002.png
    :class: sphx-glr-single-img




%% Test asymptotic behavour


.. code-block:: default

    def dip_potential(Reval, Rdip, moment):
        R  = Reval - Rdip
        r = np.linalg.norm(R, axis=1)
        return (moment*R).sum(axis=1)/r**3

    # Center of mass
    Rdip = points.mean(axis=0)
    # Moment
    m = ta[0]*tn[0]
    # Eval points
    Neval = 100
    p_eval2 = np.zeros((Neval, 3))
    z = np.linspace(0.01,100, Neval)
    p_eval2[:,2] = z
    p_eval2 += Rdip


    plt.figure()

    # Plot dipole field approximating uniform dipolar density
    plt.semilogy(z, dip_potential(p_eval2, Rdip, m))
    # Plot sum of the linear dipoles
    RR = p_eval2[:,None,None,:] - p_tris[None,:,:,:]
    pot = triangle_potential_dipole_linear(RR, tn, ta, False)
    plt.semilogy(z,  pot.sum(axis=-1)[:,0])






.. image:: /auto_examples/validation/images/sphx_glr_validate_linear_dipole_003.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  1.209 seconds)

**Estimated memory usage:**  9 MB


.. _sphx_glr_download_auto_examples_validation_validate_linear_dipole.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: validate_linear_dipole.py <validate_linear_dipole.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: validate_linear_dipole.ipynb <validate_linear_dipole.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
