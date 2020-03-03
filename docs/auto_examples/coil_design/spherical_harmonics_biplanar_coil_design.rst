.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py:


High-order spherical harmonic biplanar coil design
==================================================

Example showing a basic biplanar coil producing a high-order spherical harmonic field
in a target region between the two coil planes.



.. code-block:: default


    import numpy as np
    from mayavi import mlab
    import trimesh

    from bfieldtools.mesh_class import Conductor
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
    center_offset = np.array([0, 0, 0]) * scaling_factor
    standoff = np.array([0, 4, 0]) * scaling_factor

    #Create coil plane pairs
    coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                             planemesh.faces, process=False)

    coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                         planemesh.faces, process=False)

    joined_planes = coil_plus.union(coil_minus)

    #Create mesh class object
    coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)







Set up target and stray field points


.. code-block:: default


    #Here, the target points are on a volumetric grid within a sphere

    center = np.array([0, 0, 0]) * scaling_factor

    sidelength = 2 * scaling_factor
    n = 8
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



    #    #Here, the stray field points are on a spherical surface
    stray_radius = 20 * scaling_factor
    #    stray_length = 20 * scaling_factor
    #
    #    stray_points = cylinder_points(radius=stray_radius,
    #                                   length = stray_length,
    #                                   nlength = 5,
    #                                   nalpha = 30,
    #                                   orientation=np.array([1, 0, 0]))
    #
    stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
    stray_points = stray_points_mesh.vertices + center

    n_stray_points = len(stray_points)











Create bfield specifications used when optimizing the coil geometry


.. code-block:: default


    #The absolute target field amplitude is not of importance,
    # and it is scaled to match the C matrix in the optimization function


    from bfieldtools import sphtools


    lmax = 4
    alm = np.zeros((lmax*(lmax+2),))
    blm = np.zeros((lmax*(lmax+2),))

    #
    #alm[22]+=1
    blm[22]+=1

    sphfield = sphtools.field(target_points, alm, blm, lmax)

    target_field = sphfield/np.max(sphfield[:, 0])

    #target_field[:, 2] = 0


    coil.plot_mesh()
    mlab.quiver3d(*target_points.T, *sphfield.T)



    target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':0, 'abs_error':0.1, 'target':target_field}
    stray_spec = {'coupling':coil.B_coupling(stray_points), 'abs_error':0.01, 'rel_error':0, 'target':np.zeros((n_stray_points, 3))}

    bfield_specification = [target_spec, stray_spec]




.. image:: /auto_examples/coil_design/images/sphx_glr_spherical_harmonics_biplanar_coil_design_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 3184 vertices by 160 target points... took 0.26 seconds.
    Computing magnetic field coupling matrix, 3184 vertices by 642 target points... took 0.89 seconds.



Run QP solver


.. code-block:: default

    import mosek

    coil.j, prob = optimize_streamfunctions(coil,
                                       [target_spec, stray_spec],
                                       objective='minimum_inductive_energy',
                                       solver='MOSEK',
                                       solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                       )




.. code-block:: pytb

    Traceback (most recent call last):
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 480, in _memory_usage
        out = func()
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 465, in __call__
        exec(self.code, self.globals)
      File "/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools/examples/coil_design/spherical_harmonics_biplanar_coil_design.py", line 135, in <module>
        solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
      File "/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools/bfieldtools/coil_optimize.py", line 269, in optimize_streamfunctions
        constraints = [G@x >= lb, G@x <= ub]
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 46, in cast_op
        return binary_op(self, other)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 592, in __ge__
        return other.__le__(self)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 46, in cast_op
        return binary_op(self, other)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 581, in __le__
        return Inequality(self, other)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/constraints/nonpos.py", line 96, in __init__
        self._expr = lhs - rhs
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 46, in cast_op
        return binary_op(self, other)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 464, in __sub__
        return self + -other
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 46, in cast_op
        return binary_op(self, other)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/expressions/expression.py", line 452, in __add__
        return cvxtypes.add_expr()([self, other])
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/atoms/affine/add_expr.py", line 33, in __init__
        super(AddExpression, self).__init__(*arg_groups)
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/atoms/atom.py", line 45, in __init__
        self._shape = self.shape_from_args()
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/atoms/affine/add_expr.py", line 41, in shape_from_args
        return u.shape.sum_shapes([arg.shape for arg in self.args])
      File "/u/80/makinea1/unix/miniconda3/lib/python3.6/site-packages/cvxpy/utilities/shape.py", line 49, in sum_shapes
        len(shapes)*" %s" % tuple(shapes))
    ValueError: Cannot broadcast dimensions  (2406,) (2886,)




Plot coil windings and target points


.. code-block:: default


    N_contours = 10

    loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

    f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
    mlab.clf()

    plot_3d_current_loops(loops, colors='auto', figure=f)

    B_target = coil.B_coupling(target_points) @ coil.j

    mlab.quiver3d(*target_points.T, *B_target.T)

.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  8.318 seconds)


.. _sphx_glr_download_auto_examples_coil_design_spherical_harmonics_biplanar_coil_design.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: spherical_harmonics_biplanar_coil_design.py <spherical_harmonics_biplanar_coil_design.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: spherical_harmonics_biplanar_coil_design.ipynb <spherical_harmonics_biplanar_coil_design.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
