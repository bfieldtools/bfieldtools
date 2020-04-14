.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_validation_test_integrals.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_validation_test_integrals.py:


Integrals testing
==================================================


.. code-block:: default


    import numpy as np

    import matplotlib.pyplot as plt

    import trimesh
    from mayavi import mlab








Test potential shape slightly above the surface
########################################################


.. code-block:: default

    x = np.sin(np.pi / 6)
    y = np.cos(np.pi / 6)
    points = (
        np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [x, y, 0],
                [-x, y, 0],
                [-1, 0, 0],
                [-x, -y, 0],
                [x, -y, 0],
            ]
        )
        * 2
    )

    tris = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 6], [0, 6, 1]])
    mesh = trimesh.Trimesh(points, tris)
    scalars = np.zeros(7)
    scalars[0] = 1








Linear dipole density
########################################################


.. code-block:: default


    # Sign ok
    points = np.array([[0.1, 1, 1], [0.1, 1, -1], [0.1, -1, -1], [0.1, -1, 1]]) * 2
    # points = np.roll(points, 2, 1)
    tris = np.array([[0, 1, 2], [2, 3, 0]])
    mesh2 = trimesh.Trimesh(points, tris)
    for ii in range(7):
        mesh2 = mesh2.subdivide()

    from bfieldtools.integrals_old import triangle_potential_dipole_linear as t1
    from bfieldtools.integrals import triangle_potential_dipole_linear as t2

    RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
    p1 = t1(RR, mesh.face_normals, mesh.area_faces, planar=False)
    p2 = t2(RR, mesh.face_normals, mesh.area_faces)

    assert np.allclose(p1, p2)


    mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p1[:, :, 0].sum(axis=1))
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

    mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p2[:, :, 0].sum(axis=1))
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

    mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    mlab.triangular_mesh(
        *mesh2.vertices.T, mesh2.faces, scalars=(p1 - p2)[:, :, 0].sum(axis=1)
    )
    mlab.colorbar()




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_003.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.core.lut_manager.LUTManager object at 0x000001F8015E2E08>




.. code-block:: default

    points = np.zeros((100, 3))
    points[:, 2] = np.linspace(-1, 1, 100)
    from bfieldtools.integrals_old import omega as omega1
    from bfieldtools.integrals import omega as omega2

    RR = points[:, None, None, :] - mesh.vertices[None, mesh.faces]
    o1 = omega1(RR).sum(axis=1)
    o2 = omega2(RR).sum(axis=1)

    assert np.allclose(o1, -o2)

    plt.plot(o1)
    plt.plot(o2)
    mlab.plot3d(*points.T, points[:, 2], colormap="seismic")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)




.. image:: /auto_examples/validation/images/sphx_glr_test_integrals_004.png
    :class: sphx-glr-single-img

.. image:: /auto_examples/validation/images/sphx_glr_test_integrals_005.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x000001F80588BEB8>



Plot x_i


.. code-block:: default


    from bfieldtools.integrals import x_distance

    RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
    xdist = x_distance(RR, mesh.face_normals)
    mlab.triangular_mesh(
        *mesh2.vertices.T,
        mesh2.faces,
        scalars=xdist[:, 1, 0],
        vmin=-1,
        vmax=1,
        colormap="seismic"
    )
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)




.. image:: /auto_examples/validation/images/sphx_glr_test_integrals_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x000001F870784200>



Uniform charge density


.. code-block:: default

    from bfieldtools.integrals_old import triangle_potential_uniform as u1
    from bfieldtools.integrals import triangle_potential_uniform as u2

    RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
    p1 = u1(RR, mesh.face_normals, planar=False)
    p2 = u2(RR, mesh.face_normals, planar=False)

    assert np.allclose(p1, p2)


    mlab.figure("uniform charge density (old)")
    mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p1.sum(axis=1))
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)
    mlab.figure("uniform charge density (new)")
    mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p2.sum(axis=1))
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_007.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_008.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x000001F8709FCE60>




.. code-block:: default

    from bfieldtools.integrals import d_distance

    RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
    ddist = d_distance(RR, mesh.face_normals)
    mlab.figure("d distance")
    mlab.triangular_mesh(
        *mesh2.vertices.T,
        mesh2.faces,
        scalars=ddist[:, 0],
        vmin=-1,
        vmax=1,
        colormap="seismic"
    )
    mlab.colorbar()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
    mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)




.. image:: /auto_examples/validation/images/sphx_glr_test_integrals_009.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <mayavi.modules.vectors.Vectors object at 0x000001F80082C2B0>




.. code-block:: default

    from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic_old
    from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic

    b1 = magnetic_field_coupling_analytic_old(mesh, mesh2.vertices)
    b2 = magnetic_field_coupling_analytic(mesh, mesh2.vertices)

    assert np.allclose(b1, b2)

    mlab.figure("b field")
    mlab.quiver3d(*mesh2.vertices.T, *b1[:, :, 0].T)
    mlab.quiver3d(*mesh2.vertices.T, *b2[:, :, 0].T)




.. image:: /auto_examples/validation/images/sphx_glr_test_integrals_010.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix analytically, 7 vertices by 16641 target points... took 0.13 seconds.
    Computing magnetic field coupling matrix analytically, 7 vertices by 16641 target points... took 0.12 seconds.

    <mayavi.modules.vectors.Vectors object at 0x000001F8015CD468>



Gammma


.. code-block:: default

    from bfieldtools.integrals_old import gamma0 as g1
    from bfieldtools.integrals import gamma0 as g2

    # RR =  mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
    t = np.linspace(-1.5, 1.5)
    points = (
        t[:, None] * mesh.vertices[mesh.faces][0][0]
        + (1 - t)[:, None] * mesh.vertices[mesh.faces][0][1]
    )


    R = points[:, None, None, :] - mesh.vertices[None, mesh.faces]
    p1 = g1(R, symmetrize=True)
    p2 = g2(R, symmetrize=True)

    assert np.allclose(p1, p2)

    plt.figure()
    plt.plot(p1[:, 0, :])
    plt.figure()
    plt.plot(p2[:, 0, :])



.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_011.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/validation/images/sphx_glr_test_integrals_012.png
            :class: sphx-glr-multi-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    [<matplotlib.lines.Line2D object at 0x000001F8015AA6D8>, <matplotlib.lines.Line2D object at 0x000001F8015AA780>, <matplotlib.lines.Line2D object at 0x000001F8015AAE10>]




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  10.906 seconds)


.. _sphx_glr_download_auto_examples_validation_test_integrals.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: test_integrals.py <test_integrals.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: test_integrals.ipynb <test_integrals.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
