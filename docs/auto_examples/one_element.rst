.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_one_element.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_one_element.py:


Created on Tue Sep  3 15:46:47 2019

@author: makinea1


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
    from bfieldtools.laplacian_mesh import gradient
    from bfieldtools.magnetic_field_mesh import compute_U, compute_A

    import trimesh
    from mayavi import mlab







%% Test potential shape slightly above the surface


.. code-block:: default

    x = np.sin(np.pi/6)
    y = np.cos(np.pi/6)
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [x, y, 0],
                       [-x, y, 0],
                       [-1, 0, 0],
                       [-x, -y, 0],
                       [x, -y, 0]])

    tris = np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,1]])
    mesh = trimesh.Trimesh(points, tris)
    scalars = np.zeros(7)
    scalars[0] = 1
    s1 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')

    s2 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')
    s2.enable_contours = True
    s2.actor.mapper.scalar_range = np.array([0., 1.])
    s2.actor.mapper.scalar_visibility = False
    s2.actor.property.render_lines_as_tubes = True

    points = np.array([[0.01, 1, 1],
                       [0.01, 1, -1],
                       [0.01, -1, -1],
                       [0.01, -1, 1]])*2
    tris=np.array([[0,1,2], [2,3,0]])
    mesh2 = trimesh.Trimesh(points, tris)
    for ii in range(7):
        mesh2 =mesh2.subdivide()

    U = compute_U(mesh, mesh2.vertices) @ scalars

    s3= mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=U, colormap='bwr')
    s3.enable_contours = True
    s3.contour.minimum_contour = -5.2e-07
    s3.contour.maximum_contour = 5.2e-07
    s3.actor.property.render_lines_as_tubes = True


    points = np.array([[1, 1, 0.1],
                       [1, -1, 0.1],
                       [-1, -1, 0.1],
                       [-1, 1, 0.1]])*2
    tris=np.array([[0,1,2], [2,3,0]])
    mesh3 = trimesh.Trimesh(points, tris)
    for ii in range(7):
        mesh3 =mesh3.subdivide()
    A = compute_A(mesh, mesh3.vertices) @ scalars
    mlab.quiver3d(*mesh3.vertices.T, *A)


.. image:: /auto_examples/images/sphx_glr_one_element_001.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing potential matrix
    0 % computed
    1 % computed
    1 % computed
    2 % computed
    3 % computed
    3 % computed
    4 % computed
    4 % computed
    5 % computed
    6 % computed
    6 % computed
    7 % computed
    7 % computed
    8 % computed
    9 % computed
    9 % computed
    10 % computed
    10 % computed
    11 % computed
    12 % computed
    12 % computed
    13 % computed
    13 % computed
    14 % computed
    15 % computed
    15 % computed
    16 % computed
    16 % computed
    17 % computed
    18 % computed
    18 % computed
    19 % computed
    20 % computed
    20 % computed
    21 % computed
    21 % computed
    22 % computed
    23 % computed
    23 % computed
    24 % computed
    24 % computed
    25 % computed
    26 % computed
    26 % computed
    27 % computed
    27 % computed
    28 % computed
    29 % computed
    29 % computed
    30 % computed
    30 % computed
    31 % computed
    32 % computed
    32 % computed
    33 % computed
    33 % computed
    34 % computed
    35 % computed
    35 % computed
    36 % computed
    36 % computed
    37 % computed
    38 % computed
    38 % computed
    39 % computed
    39 % computed
    40 % computed
    41 % computed
    41 % computed
    42 % computed
    42 % computed
    43 % computed
    44 % computed
    44 % computed
    45 % computed
    45 % computed
    46 % computed
    47 % computed
    47 % computed
    48 % computed
    48 % computed
    49 % computed
    50 % computed
    50 % computed
    51 % computed
    51 % computed
    52 % computed
    53 % computed
    53 % computed
    54 % computed
    54 % computed
    55 % computed
    56 % computed
    56 % computed
    57 % computed
    57 % computed
    58 % computed
    59 % computed
    59 % computed
    60 % computed
    60 % computed
    61 % computed
    62 % computed
    62 % computed
    63 % computed
    63 % computed
    64 % computed
    65 % computed
    65 % computed
    66 % computed
    66 % computed
    67 % computed
    68 % computed
    68 % computed
    69 % computed
    69 % computed
    70 % computed
    71 % computed
    71 % computed
    72 % computed
    72 % computed
    73 % computed
    74 % computed
    74 % computed
    75 % computed
    75 % computed
    76 % computed
    77 % computed
    77 % computed
    78 % computed
    78 % computed
    79 % computed
    80 % computed
    80 % computed
    81 % computed
    81 % computed
    82 % computed
    83 % computed
    83 % computed
    84 % computed
    84 % computed
    85 % computed
    86 % computed
    86 % computed
    87 % computed
    87 % computed
    88 % computed
    89 % computed
    89 % computed
    90 % computed
    90 % computed
    91 % computed
    92 % computed
    92 % computed
    93 % computed
    93 % computed
    94 % computed
    95 % computed
    95 % computed
    96 % computed
    96 % computed
    97 % computed
    98 % computed
    98 % computed
    99 % computed
    100 % computed
    Computing potential matrix
    0 % computed
    1 % computed
    1 % computed
    2 % computed
    3 % computed
    3 % computed
    4 % computed
    4 % computed
    5 % computed
    6 % computed
    6 % computed
    7 % computed
    7 % computed
    8 % computed
    9 % computed
    9 % computed
    10 % computed
    10 % computed
    11 % computed
    12 % computed
    12 % computed
    13 % computed
    13 % computed
    14 % computed
    15 % computed
    15 % computed
    16 % computed
    16 % computed
    17 % computed
    18 % computed
    18 % computed
    19 % computed
    20 % computed
    20 % computed
    21 % computed
    21 % computed
    22 % computed
    23 % computed
    23 % computed
    24 % computed
    24 % computed
    25 % computed
    26 % computed
    26 % computed
    27 % computed
    27 % computed
    28 % computed
    29 % computed
    29 % computed
    30 % computed
    30 % computed
    31 % computed
    32 % computed
    32 % computed
    33 % computed
    33 % computed
    34 % computed
    35 % computed
    35 % computed
    36 % computed
    36 % computed
    37 % computed
    38 % computed
    38 % computed
    39 % computed
    39 % computed
    40 % computed
    41 % computed
    41 % computed
    42 % computed
    42 % computed
    43 % computed
    44 % computed
    44 % computed
    45 % computed
    45 % computed
    46 % computed
    47 % computed
    47 % computed
    48 % computed
    48 % computed
    49 % computed
    50 % computed
    50 % computed
    51 % computed
    51 % computed
    52 % computed
    53 % computed
    53 % computed
    54 % computed
    54 % computed
    55 % computed
    56 % computed
    56 % computed
    57 % computed
    57 % computed
    58 % computed
    59 % computed
    59 % computed
    60 % computed
    60 % computed
    61 % computed
    62 % computed
    62 % computed
    63 % computed
    63 % computed
    64 % computed
    65 % computed
    65 % computed
    66 % computed
    66 % computed
    67 % computed
    68 % computed
    68 % computed
    69 % computed
    69 % computed
    70 % computed
    71 % computed
    71 % computed
    72 % computed
    72 % computed
    73 % computed
    74 % computed
    74 % computed
    75 % computed
    75 % computed
    76 % computed
    77 % computed
    77 % computed
    78 % computed
    78 % computed
    79 % computed
    80 % computed
    80 % computed
    81 % computed
    81 % computed
    82 % computed
    83 % computed
    83 % computed
    84 % computed
    84 % computed
    85 % computed
    86 % computed
    86 % computed
    87 % computed
    87 % computed
    88 % computed
    89 % computed
    89 % computed
    90 % computed
    90 % computed
    91 % computed
    92 % computed
    92 % computed
    93 % computed
    93 % computed
    94 % computed
    95 % computed
    95 % computed
    96 % computed
    96 % computed
    97 % computed
    98 % computed
    98 % computed
    99 % computed
    100 % computed




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  2.986 seconds)

**Estimated memory usage:**  138 MB


.. _sphx_glr_download_auto_examples_one_element.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: one_element.py <one_element.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: one_element.ipynb <one_element.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
