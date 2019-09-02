.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_shielded_cylindrical_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_shielded_cylindrical_coil_design.py:


Shielded cylindrical coil
=========================
Compact example of design of a shielded cylindrical coil, for which the transient field
due to inductive interaction with the shield is minimized


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.magnetic_field_mesh import compute_C
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load example coil mesh that is centered on the origin
    coilmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True)

    coilmesh.apply_scale(0.75)

    coilmesh.vertices, coilmesh.faces = trimesh.remesh.subdivide(coilmesh.vertices, coilmesh.faces)

    #Specify offset from origin
    center_offset = np.array([0, 0, 0.75])

    #Apply offset
    coilmesh = trimesh.Trimesh(coilmesh.vertices + center_offset,
                                coilmesh.faces, process=False)

    #Create mesh class object
    coil = MeshWrapper(verts=coilmesh.vertices, tris=coilmesh.faces, fix_normals=True)

    # Separate object for shield geometry
    shield = MeshWrapper(mesh_file=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True, fix_normals=True)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!
    face_normals didn't match triangles, ignoring!



Set up target  points and plot geometry


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 3])

    sidelength = 0.75 * scaling_factor
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


    #Plot coil, shield and target points

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                    size=(800, 800))

    coil.plot_mesh()
    shield.plot_mesh()
    mlab.points3d(*target_points.T)







.. image:: /auto_examples/coil_design/images/sphx_glr_shielded_cylindrical_coil_design_001.png
    :class: sphx-glr-single-img




Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)
    shield.C = compute_C(shield.mesh, target_points)

    mutual_inductance = mutual_inductance_matrix(coil.mesh, shield.mesh)

    # Take into account the field produced by currents induced into the shield
    # NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

    shield.coupling = -np.linalg.pinv(shield.inductance) @ mutual_inductance.T
    secondary_C = (shield.C.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3536 vertices by 672 target points... took 1.02 seconds.
    Computing C matrix, 904 vertices by 672 target points... took 0.21 seconds.
    Calculating potentials
    Inserting stuff into M-matrix
    Computing inductance matrix in 1 chunks since 9 GiB memory is available...
    Calculating potentials, chunk 1/1
    Inductance matrix computation took 6.06 seconds.



Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1 # Homogeneous Z-field

    target_spec = {'C':coil.C, 'rel_error':0.01, 'abs_error':0, 'target_field':target_field}


    induction_spec = {'C':secondary_C, 'abs_error':0.1, 'rel_error':0, 'target_field':np.zeros(target_field.shape)}







Run QP solver


.. code-block:: default


    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 0.5

    coil.I, coil.sol = optimize_streamfunctions(coil,
                                                [target_spec, induction_spec],
                                                laplacian_smooth=0,
                                                tolerance=tolerance)

    shield.induced_I = shield.coupling @ coil.I






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 9 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 96.50 seconds.
    Scaling matrices before optimization. This requires singular value computation, hold on.
    Solving quadratic programming problem using cvxopt...
         pcost       dcost       gap    pres   dres
     0:  4.5366e+01  7.7593e+02  1e+04  3e+00  4e-14
     1:  1.0096e+02  1.1497e+03  3e+03  9e-01  3e-14
     2:  5.6082e+02  2.4779e+03  3e+03  6e-01  9e-14
     3:  5.7943e+02  2.6428e+03  3e+03  6e-01  9e-14
     4:  1.0195e+03  5.8843e+03  3e+03  5e-01  2e-13
     5:  2.6725e+03  1.0689e+04  4e+03  4e-01  4e-13
    Optimal solution found.



Plot coil windings and target points


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    surface = mlab.pipeline.triangular_mesh_source(*coil.mesh.vertices.T, coil.mesh.faces,scalars=coil.I)

    windings = mlab.pipeline.contour_surface(surface, contours=20)
    windings.module_manager.scalar_lut_manager.number_of_colors = 2 #Color windings according to current direction
    windings.module_manager.scalar_lut_manager.reverse_lut = True #Flip LUT for the colors to correspond to RdBu colormap

    shield_surface = mlab.pipeline.triangular_mesh_source(*shield.mesh.vertices.T, shield.mesh.faces,scalars=shield.induced_I)

    shield_surface_render = mlab.pipeline.surface(shield_surface, colormap='RdBu')

    shield_surface_render.actor.property.frontface_culling = True

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)

    mlab.title('Coils which minimize the transient effects of conductive shield')





.. image:: /auto_examples/coil_design/images/sphx_glr_shielded_cylindrical_coil_design_002.png
    :class: sphx-glr-single-img




For comparison, let's see how the coils look when we ignore the conducting shield


.. code-block:: default



    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 0.5

    coil.unshielded_I, coil.unshielded_sol = optimize_streamfunctions(coil,
                                                [target_spec],
                                                laplacian_smooth=0,
                                                tolerance=tolerance)

    shield.unshielded_induced_I = shield.coupling @ coil.unshielded_I

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                    size=(800, 800))
    mlab.clf()

    surface = mlab.pipeline.triangular_mesh_source(*coil.mesh.vertices.T, coil.mesh.faces,scalars=coil.unshielded_I)

    windings = mlab.pipeline.contour_surface(surface, contours=20)
    windings.module_manager.scalar_lut_manager.number_of_colors = 2 #Color windings according to current direction
    windings.module_manager.scalar_lut_manager.reverse_lut = True #Flip LUT for the colors to correspond to RdBu colormap

    shield_surface = mlab.pipeline.triangular_mesh_source(*shield.mesh.vertices.T, shield.mesh.faces,scalars=shield.unshielded_induced_I)

    shield_surface_render = mlab.pipeline.surface(shield_surface, colormap='RdBu')

    shield_surface_render.actor.property.frontface_culling = True

    B_target = coil.C.transpose([0, 2, 1]) @ coil.unshielded_I

    mlab.quiver3d(*target_points.T, *B_target.T)

    mlab.title('Coils which ignore the conductive shield')




.. image:: /auto_examples/coil_design/images/sphx_glr_shielded_cylindrical_coil_design_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Scaling matrices before optimization. This requires singular value computation, hold on.
    Solving quadratic programming problem using cvxopt...
         pcost       dcost       gap    pres   dres
     0:  4.3937e+01  6.6183e+01  5e+03  2e+00  4e-14
     1:  5.0479e+01  6.8140e+01  2e+02  1e-01  2e-14
     2:  6.1947e+01  9.9187e+01  1e+02  4e-02  3e-14
     3:  6.7204e+01  1.2054e+02  1e+02  3e-02  5e-14
     4:  7.3748e+01  1.8486e+02  7e+01  2e-02  1e-13
    Optimal solution found.




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  53.747 seconds)


.. _sphx_glr_download_auto_examples_coil_design_shielded_cylindrical_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: shielded_cylindrical_coil_design.py <shielded_cylindrical_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: shielded_cylindrical_coil_design.ipynb <shielded_cylindrical_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
