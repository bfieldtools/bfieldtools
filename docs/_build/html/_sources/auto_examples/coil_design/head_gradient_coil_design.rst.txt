.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_head_gradient_coil_design.py:


Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet


.. code-block:: default



    import numpy as np
    import matplotlib.pyplot as plt
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    helmetmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools',
                                                                       'example_meshes/meg_helmet.obj'),
                              process=False)

    #planemesh.apply_scale(scaling_factor)
    #
    ##Specify coil plane geometry
    #center_offset = np.array([0, 0, 0]) * scaling_factor
    #standoff = np.array([0, 4, 0]) * scaling_factor
    #
    ##Create coil plane pairs
    #coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
    #                         planemesh.faces, process=False)
    ##
    #coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
    #                     planemesh.faces, process=False)

    #joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    shapely.geometry.Polygon not available!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/creation.py", line 22, in <module>
        from shapely.geometry import Polygon
    ModuleNotFoundError: No module named 'shapely'



Set up target and stray field points.
Here, the target points are on a volumetric grid within a sphere


.. code-block:: default


    center = np.array([0, 0, 0.04]) * scaling_factor

    sidelength = 0.1 * scaling_factor
    n = 12
    xx = np.linspace(-sidelength/2, sidelength/2, n)
    yy = np.linspace(-sidelength/2, sidelength/2, n)
    zz = np.linspace(-sidelength/2, sidelength/2, n)
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    target_points = np.array([x, y, z]).T

    #Turn cube into sphere by rejecting points "in the corners"
    target_points = target_points[np.linalg.norm(target_points, axis=1) < sidelength/2]  + center








Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 2044 vertices by 672 target points... took 0.61 seconds.



Specify target field and run solver


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.ones(target_points.shape[0], ) * target_points[:,0]



    target_field = np.zeros(target_points.shape)
    target_field[:, 0] = target_field[:, 0] + 1 * target_points[:,0]/np.max(target_points[:,0])

    target_spec = {'C':coil.C, 'rel_error':0.1, 'abs_error':0, 'target_field':target_field}

    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 5

    coil.I, coil.sol = optimize_streamfunctions(coil,
                                                [target_spec],
                                                laplacian_smooth=0,
                                                tolerance=tolerance)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 1 chunks since 10 GiB memory is available...
    Calculating potentials, chunk 1/1
    Inductance matrix computation took 32.43 seconds.
    Scaling matrices before optimization. This requires singular value computation, hold on.
    Solving quadratic programming problem using cvxopt...
         pcost       dcost       gap    pres   dres
     0:  2.9677e+01  3.3517e+02  1e+04  5e+00  2e-14
     1:  6.0433e+01  6.1774e+02  4e+03  2e+00  3e-14
     2:  1.0599e+02  1.5067e+03  3e+03  1e+00  5e-14
    Optimal solution found.



Plot coil windings and magnetic field in target points


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    surface = mlab.pipeline.triangular_mesh_source(*coil.mesh.vertices.T, coil.mesh.faces,scalars=coil.I)

    windings = mlab.pipeline.contour_surface(surface, contours=10)


    B_target = np.vstack((coil.C[:, :, 0].dot(coil.I),
                      coil.C[:, :, 1].dot(coil.I),
                      coil.C[:, :, 2].dot(coil.I))).T


    mlab.quiver3d(*target_points.T, *B_target.T)
    f.scene.isometric_view()



.. image:: /auto_examples/coil_design/images/sphx_glr_head_gradient_coil_design_001.png
    :class: sphx-glr-single-img





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  38.254 seconds)


.. _sphx_glr_download_auto_examples_coil_design_head_gradient_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: head_gradient_coil_design.py <head_gradient_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: head_gradient_coil_design.ipynb <head_gradient_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
