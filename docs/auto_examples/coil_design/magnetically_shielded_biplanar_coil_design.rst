.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py:


Magnetically shielded  coil
===========================
Compact example of design of a biplanar coil within a cylindrical shield.
The effect of the shield is prospectively taken into account while designing the coil.
The coil is positioned close to the end of the shield to demonstrate the effect


.. code-block:: default



    import numpy as np
    from mayavi import mlab
    import trimesh


    from bfieldtools.mesh_class import MeshWrapper
    from bfieldtools.magnetic_field_mesh import compute_C, compute_U
    from bfieldtools.coil_optimize import optimize_streamfunctions
    from bfieldtools.contour import scalar_contour
    from bfieldtools.viz import plot_3d_current_loops

    import pkg_resources


    #Set unit, e.g. meter or millimeter.
    # This doesn't matter, the problem is scale-invariant
    scaling_factor = 1


    #Load simple plane mesh that is centered on the origin
    planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

    planemesh.apply_scale(scaling_factor)

    #Specify coil plane geometry
    center_offset = np.array([9, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = MeshWrapper(mesh_obj=joined_planes, fix_normals=True)

    # Separate object for shield geometry
    shieldmesh = trimesh.load('/l/bfieldtools/bfieldtools/example_meshes/closed_cylinder.stl')
    shieldmesh.apply_scale(15)

    shield = MeshWrapper(mesh_obj=shieldmesh, process=True, fix_normals=True)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!



Set up target  points and plot geometry


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere
    # Set up target and stray field points

    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([9, 0, 0]) * scaling_factor

    sidelength = 3 * scaling_factor
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

    coil.plot_mesh(representation='surface')
    shield.plot_mesh()
    mlab.points3d(*target_points.T)





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_001.png
    :class: sphx-glr-single-img




Compute C matrices that are used to compute the generated magnetic field


.. code-block:: default


    coil.C = compute_C(coil.mesh, target_points)
    shield.C = compute_C(shield.mesh, target_points)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing C matrix, 3184 vertices by 672 target points... took 0.97 seconds.
    Computing C matrix, 962 vertices by 672 target points... took 0.27 seconds.



Let's design a coil without taking the magnetic shield into account


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function
    target_field = np.zeros(target_points.shape)
    target_field[:, 1] = target_field[:, 1] + 1 # Homogeneous Z-field

    target_spec = {'C':coil.C, 'rel_error':0.01, 'abs_error':0, 'target_field':target_field}


    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 0.5

    coil.I, coil.sol = optimize_streamfunctions(coil,
                                                [target_spec],
                                                objective='minimum_inductive_energy',
                                                tolerance=tolerance)








.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing inductance matrix in 2 chunks since 7 GiB memory is available...
    Calculating potentials, chunk 1/2
    Calculating potentials, chunk 2/2
    Inductance matrix computation took 67.63 seconds.
    Solving quadratic programming problem using cvxopt...
         pcost       dcost       gap    pres   dres
     0:  4.3858e+01  7.6264e+01  5e+03  2e+00  3e-14
     1:  5.2944e+01  9.1403e+01  6e+02  2e-01  2e-14
     2:  1.0614e+02  1.3543e+02  2e+02  6e-02  3e-14
     3:  1.0230e+02  1.7052e+02  2e+02  5e-02  5e-14
     4:  1.0406e+02  1.9056e+02  2e+02  5e-02  9e-14
     5:  1.8354e+02  3.6462e+02  2e+02  3e-02  5e-13
     6:  1.9672e+02  4.7650e+02  2e+02  2e-02  5e-13
    Optimal solution found.



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=10)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_002.png
    :class: sphx-glr-single-img




Now, let's compute the effect of the shield on the field produced by the coil


.. code-block:: default


    # Calculate primary potential matrix at the shield surface
    P_prim = compute_U(coil.mesh, shield.mesh.vertices)

    # Calculate linear collocation BEM matrix
    P_bem = compute_U(shield.mesh, shield.mesh.vertices)

    # Recalculate diag elements according to de Munck paper
    for diag_index in range(P_bem.shape[0]):
        P_bem[diag_index, diag_index] = 0
        P_bem[diag_index, diag_index] = -P_bem[diag_index, :].sum()

    # Matrix misses one rank, make it invertible
    # by rank-one update (sets potential of constant dipole layer)
    P_bem += np.ones(P_bem.shape)/P_bem.shape[0]


    # Solve equivalent stream function for the perfect linear mu-metal layer.
    # This is the equivalent surface current in the shield that would cause its
    # scalar magnetic potential to be constant
    shield.I =  np.linalg.solve(P_bem, P_prim @ coil.I)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing U matrix, 3184 vertices by 962 target points... took 12.59 seconds.
    Computing U matrix, 962 vertices by 962 target points... took 3.97 seconds.



Plot the difference in field when taking the shield into account


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    B_target_w_shield = coil.C.transpose([0, 2, 1]) @ coil.I + shield.C.transpose([0, 2, 1]) @ shield.I

    B_quiver = mlab.quiver3d(*target_points.T, *(B_target_w_shield - B_target).T, colormap='viridis', mode='arrow')
    f.scene.isometric_view()
    mlab.colorbar(B_quiver, title='Difference in magnetic field (a.u.)')

    import seaborn as sns
    import matplotlib.pyplot as plt




    fig, axes = plt.subplots(1, 3, figsize=(10, 4))

    fig.suptitle('Component-wise effect of magnetic shield on target field amplitude distribution')
    for ax_idx, ax in enumerate(axes):

        sns.distplot(B_target[:, ax_idx], label='Without shield', ax=ax)
        sns.distplot(B_target_w_shield[:, ax_idx], label='With shield', ax=ax)
        ax.set_xlabel('Magnetic field (a.u.)')

        if ax_idx == 2:
            ax.legend()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])





.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_003.png
    :class: sphx-glr-single-img

.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_004.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    This object has no scalar data



Let's redesign the coil taking the shield into account prospectively


.. code-block:: default


    shield.coupling = np.linalg.solve(P_bem, P_prim)

    secondary_C = (shield.C.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))

    total_C = coil.C + secondary_C

    target_spec_w_shield = {'C':total_C, 'rel_error':0.01, 'abs_error':0, 'target_field':target_field}


    # The tolerance parameter will determine the spatial detail of the coil.
    # Smaller tolerance means better but more intricate patterns. Too small values
    # will not be solveable.
    tolerance = 0.5

    coil.I2, coil.sol2 = optimize_streamfunctions(coil,
                                                [target_spec_w_shield],
                                                objective='minimum_inductive_energy',
                                                tolerance=tolerance)






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Solving quadratic programming problem using cvxopt...
         pcost       dcost       gap    pres   dres
     0:  4.1091e+01  6.7938e+01  5e+03  2e+00  4e-14
     1:  4.7921e+01  7.7016e+01  5e+02  2e-01  2e-14
     2:  8.4260e+01  1.1282e+02  2e+02  5e-02  3e-14
     3:  8.2330e+01  1.4557e+02  2e+02  4e-02  4e-14
     4:  8.3642e+01  1.6328e+02  2e+02  4e-02  9e-14
     5:  1.5115e+02  3.4381e+02  2e+02  2e-02  5e-13
     6:  1.4870e+02  3.7736e+02  2e+02  2e-02  5e-13
    Optimal solution found.



Plot coil windings and target points


.. code-block:: default


    loops, loop_values= scalar_contour(coil.mesh, coil.I2, N_contours=10)
    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.C.transpose([0, 2, 1]) @ coil.I

    mlab.quiver3d(*target_points.T, *B_target.T)

    B_target2 = coil.C.transpose([0, 2, 1]) @ coil.I2


    mlab.quiver3d(*target_points.T, *B_target2.T)




.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_005.png
    :class: sphx-glr-single-img




Finally, plot the difference in stream functions


.. code-block:: default


    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    RE_I = mlab.triangular_mesh(*coil.mesh.vertices.T, coil.mesh.faces, scalars=100 * (coil.I-coil.I2)/coil.I, colormap='RdBu')
    mlab.colorbar(RE_I, title='Relative error (%)')


.. image:: /auto_examples/coil_design/images/sphx_glr_magnetically_shielded_biplanar_coil_design_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    /l/bfieldtools/examples/coil_design/magnetically_shielded_biplanar_coil_design.py:249: RuntimeWarning: invalid value encountered in true_divide
      RE_I = mlab.triangular_mesh(*coil.mesh.vertices.T, coil.mesh.faces, scalars=100 * (coil.I-coil.I2)/coil.I, colormap='RdBu')




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 2 minutes  5.402 seconds)

**Estimated memory usage:**  8033 MB


.. _sphx_glr_download_auto_examples_coil_design_magnetically_shielded_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: magnetically_shielded_biplanar_coil_design.py <magnetically_shielded_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: magnetically_shielded_biplanar_coil_design.ipynb <magnetically_shielded_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
