.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_laplacian_w_holes.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_laplacian_w_holes.py:


Created on Mon Nov  4 09:37:19 2019

@author: Antti


.. code-block:: default


    import numpy as np
    import trimesh
    from timeit import timeit
    from mayavi import mlab
    import matplotlib.pyplot as plt

    import pkg_resources
    import sys
    #path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
    #if path not in sys.path:
    #    sys.path.insert(0,path)

    from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix
    from bfieldtools.mesh_class import MeshWrapper


    #Load simple plane mesh that is centered on the origin
    file_obj = pkg_resources.resource_filename('bfieldtools',
                        'example_meshes/plane_w_holes.stl')
    coilmesh = trimesh.load(file_obj, process=True)
    coil = MeshWrapper(mesh_obj = coilmesh)

    mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces)
    mlab.points3d([-6,6,0], [0,0,0], [0,0,0])
    #mlab.show()

    #%%
    i1 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([6,0,0]))**2, axis=1) < 3**2
    i2 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([0,0,0]))**2, axis=1) < 3**2
    i3 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([-6,0,0]))**2, axis=1) < 3**2

    b1 = coil.boundary_verts[i1]
    b2 = coil.boundary_verts[i2]
    b3 = coil.boundary_verts[i3]
    b4 = coil.boundary_verts[np.logical_not(i1+i2+i3)]

    #%%
    L = laplacian_matrix(coilmesh)
    M = mass_matrix(coilmesh)
    Linner = L[coil.inner_verts, :][:, coil.inner_verts]
    Minner = M[coil.inner_verts, :][:, coil.inner_verts]

    L1 = np.array(np.sum(L[b1,:][:, coil.inner_verts], axis=0))
    L2 = np.array(np.sum(L[b2,:][:, coil.inner_verts], axis=0))
    L3 = np.array(np.sum(L[b3,:][:, coil.inner_verts], axis=0))

    u1 = np.linalg.solve(Linner.toarray() , -L1[0])

    scalars = np.zeros(L.shape[0])
    scalars[coil.inner_verts] =u1
    mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=scalars)
    #mlab.show()


    #%%
    L_holes = np.concatenate((Linner.toarray(), L1.T, L2.T, L3.T), axis=1)
    L11 = np.concatenate((L1, np.array([[-L1.sum(), 0, 0]])), axis=1)
    L21 = np.concatenate((L2, np.array([[0, -L2.sum(), 0]])), axis=1)
    L31 = np.concatenate((L3, np.array([[0, 0, -L3.sum()]])), axis=1)
    L_holes = np.concatenate((L_holes, L11, L21, L31) ,axis=0)

    #%%
    m = M.diagonal()
    M_holes = np.diag(np.concatenate((Minner.diagonal(), np.array([np.sum(m[b1]),
                                     np.sum(m[b2]), np.sum(m[b3])]))))

    #%%
    from scipy.linalg import eigh
    ee, vv = eigh(-L_holes, M_holes, eigvals=(0,10))




.. code-block:: pytb

    Traceback (most recent call last):
      File "/l/conda-envs/mne/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 474, in _memory_usage
        multiprocess=True)
      File "/l/conda-envs/mne/lib/python3.6/site-packages/memory_profiler.py", line 336, in memory_usage
        returned = f(*args, **kw)
      File "/l/conda-envs/mne/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 465, in __call__
        exec(self.code, self.globals)
      File "/l/bfieldtools/examples/laplacian_w_holes.py", line 20, in <module>
        from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix
      File "/l/bfieldtools/bfieldtools/__init__.py", line 9, in <module>
        from . import thermal_noise
      File "/l/bfieldtools/bfieldtools/thermal_noise.py", line 11, in <module>
        from .laplacian_mesh import laplacian_matrix, mass_matrix, laplacian_matrix_w_holes, mass_matrix_w_holes
    ModuleNotFoundError: No module named 'bfieldtools.laplacian_mesh'




%%


.. code-block:: default


    from scipy.linalg import eigh
    u, v = eigh(-L_holes, M_holes)


    #Normalize the laplacien eigenvectors

    for i in range(v.shape[1]):
        v[:, i] = v[:, i]/np.sqrt(u[i])

    #Assign values per vertex
    vl = np.zeros((M.shape[0], v.shape[1]))

    vl[coil.inner_verts] = v[:-3]
    vl[b1] = v[-3]
    vl[b2] = v[-2]
    vl[b3] = v[-1]


    verts = coil.mesh.vertices
    tris = coil.mesh.faces

    Nmodes = 600
    scale = 1
    contours=True
    colormap='bwr'

    mlab.figure()

    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)

        x = scale*verts[:, 0] + i*(np.max(verts[:, 0]) - np.min(verts[:, 0]))*1.2
        y = scale*verts[:, 1]+ j*(np.max(verts[:, 1]) - np.min(verts[:, 1]))*1.2
        z = scale*verts[:, 2]

        limit = np.max(np.abs(vl[:, ii]))

        s = mlab.triangular_mesh(x, y, z, tris, scalars=vl[:, ii], colormap=colormap)

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = contours


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  0.215 seconds)

**Estimated memory usage:**  9 MB


.. _sphx_glr_download_auto_examples_laplacian_w_holes.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: laplacian_w_holes.py <laplacian_w_holes.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: laplacian_w_holes.ipynb <laplacian_w_holes.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
