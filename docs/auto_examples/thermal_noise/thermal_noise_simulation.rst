.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_thermal_noise_thermal_noise_simulation.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_thermal_noise_thermal_noise_simulation.py:


Thermal noise computation
==========================

Three different examples:
   unit_sphere: DC Bnoise of a spherical shell at origin and comparison to analytical formula
   unit_disc: DC Bnoise of a unit disc at z-axis and comparison to analytical formula
   AC: AC Bnoise of a unit disc at one position


.. code-block:: default



    import numpy as np
    import matplotlib.pyplot as plt
    import trimesh
    from mayavi import mlab

    from bfieldtools.mesh_properties import self_inductance_matrix, resistance_matrix
    from bfieldtools.thermal_noise import (
        compute_current_modes_ind_res,
        noise_covar,
        noise_var,
        visualize_current_modes,
    )
    from bfieldtools.mesh_magnetics import magnetic_field_coupling

    import pkg_resources


    font = {"family": "normal", "weight": "normal", "size": 16}
    plt.rc("font", **font)

    # Fix the simulation parameters
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4 * np.pi * 1e-7
    freqs = np.array((0,))


    Nchunks = 8
    quad_degree = 2








Unit sphere
------------


.. code-block:: default



    Np = 10
    radius = np.linspace(0.1, 1, Np)
    fp = np.zeros((1, 3))

    B = np.zeros((Np, 3))
    for i in range(Np):
        mesh = trimesh.load(
            pkg_resources.resource_filename("bfieldtools", "example_meshes/unit_sphere.stl")
        )
        mesh.apply_scale(radius[i])

        B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

        S = np.ones(mesh.triangles_center.shape[0]) * sigma
        sheet_resistance = 1 / (d * S)

        # Compute the resistance and inductance matrices
        R = resistance_matrix(mesh, sheet_resistance=sheet_resistance)
        M = self_inductance_matrix(mesh, Nchunks=Nchunks, quad_degree=quad_degree)

        vl = compute_current_modes_ind_res(mesh, M, R, freqs, T, closed=True)

        #    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
        #               size=(800, 800))
        #    visualize_current_modes(mesh,vl[:,:,0], 8, 1)

        #    vl[:,0] = np.zeros(vl[:,0].shape) # fix DC-component

        Btemp = noise_var(mesh, B_coupling, vl)
        #    Btemp = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)
        B[i] = Btemp[:, :, 0]

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
    scene.scene.z_minus_view()
    surface = scene.children[0].children[0].children[0].children[0]
    surface.actor.property.representation = "wireframe"
    surface.actor.mapper.scalar_visibility = False
    scene.scene.camera.position = [0.0, 0.0, -5.530686305704514]
    scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.0, 1.0, 0.0]
    scene.scene.camera.clipping_range = [3.485379442647469, 8.118646600290083]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    scene.scene.camera.position = [0.0, 0.0, -4.570815128681416]
    scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.0, 1.0, 0.0]
    scene.scene.camera.clipping_range = [2.535106977394602, 7.1443773556116374]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    mlab.savefig(
        "/Users/joonas/Documents/Manuscripts/ThermalNoise/figures/validation/sphere.png",
        size=(800, 800),
    )

    Ban = mu0 * np.sqrt(2 * sigma * d * kB * T / (3 * np.pi * (radius) ** 2))

    plt.figure(figsize=(5, 5))
    plt.semilogy(radius, Ban * 1e15, linewidth=2, label="Analytic")
    plt.semilogy(
        radius,
        np.sqrt(B[:, 2]) * 1e15,
        "x",
        markersize=10,
        markeredgewidth=2,
        label="Numerical",
    )
    plt.grid()
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    plt.legend(frameon=False)
    plt.xlabel("Sphere radius")
    plt.ylabel(r"$B_z$ noise at DC (fT/rHz)")
    plt.tight_layout()


    RE = np.abs((np.sqrt(B[:, 2]) - Ban)) / np.abs(Ban) * 100
    plt.figure()
    plt.plot(radius, RE)
    plt.xlabel("Sphere radius")
    plt.ylabel("Relative error (%)")




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_001.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_002.png
            :class: sphx-glr-multi-img

.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_003.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    Computing magnetic field coupling matrix analytically, 2562 vertices by 1 target points... took 0.02 seconds.
    Computing self-inductance matrix using rough quadrature (degree=2).              For higher accuracy, set quad_degree to 4 or more.
    Computing 1/r-potential matrix
    findfont: Font family ['normal'] not found. Falling back to DejaVu Sans.

    Text(0, 0.5, 'Relative error (%)')



Unit disc, DC noise
---------------------


.. code-block:: default


    mesh = trimesh.load(
        pkg_resources.resource_filename("bfieldtools", "example_meshes/unit_disc.stl")
    )
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    vl = compute_current_modes(mesh)

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

    visualize_current_modes(mesh, vl, 42, 5, contours=False)

    Np = 30

    z = np.linspace(0.1, 1, Np)
    fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

    B = compute_dc_Bnoise(mesh, vl, fp, sigma, d, T)

    r = 1
    Ban = (
        mu0
        * np.sqrt(sigma * d * kB * T / (8 * np.pi * z ** 2))
        * (1 / (1 + z ** 2 / r ** 2))
    )

    plt.figure()
    plt.semilogy(z, Ban, label="Analytic")
    plt.semilogy(z, B[:, 2], "x", label="Numerical")
    plt.legend()
    plt.xlabel("Distance d/R")
    plt.ylabel("DC noise Bz (T/rHz)")
    plt.tight_layout()

    plt.figure()
    plt.plot(z, np.abs((B[:, 2] - Ban)) / np.abs(Ban) * 100)
    plt.xlabel("Distance d/R")
    plt.ylabel("Relative error (%)")



.. rst-class:: sphx-glr-script-out


.. code-block:: pytb

    Traceback (most recent call last):
      File "D:\Anaconda3\lib\site-packages\sphinx_gallery\gen_rst.py", line 460, in _memory_usage
        out = func()
      File "D:\Anaconda3\lib\site-packages\sphinx_gallery\gen_rst.py", line 442, in __call__
        exec(self.code, self.fake_main.__dict__)
      File "C:\Users\Rasmus Zetter\Documents\Aalto\bfieldtools\examples\thermal_noise\thermal_noise_simulation.py", line 144, in <module>
        vl = compute_current_modes(mesh)
    NameError: name 'compute_current_modes' is not defined




Closed cylinder, DC noise
--------------------------


.. code-block:: default


    mesh = trimesh.load(
        pkg_resources.resource_filename("bfieldtools", "example_meshes/closed_cylinder.stl")
    )
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)


    S = np.ones(mesh.triangles_center.shape[0]) * sigma
    sheet_resistance = 1 / (d * S)

    # Compute the resistance and inductance matrices
    R = resistance_matrix(mesh, sheet_resistance=sheet_resistance)
    M = self_inductance_matrix(mesh, Nchunks=Nchunks, quad_degree=quad_degree)

    vl = compute_current_modes_ind_res(mesh, M, R, freqs, T, closed=True)

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

    visualize_current_modes(mesh, vl[:, :, 0], 8, 1)


    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
    s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
    scene.scene.z_minus_view()
    surface = scene.children[0].children[0].children[0].children[0]
    surface.actor.property.representation = "wireframe"
    surface.actor.mapper.scalar_visibility = False
    scene.scene.isometric_view()
    # scene.scene.camera.position = [2.2578932293957665, 2.2578932293957665, 2.2578932293957665]
    # scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
    # scene.scene.camera.view_angle = 30.0
    # scene.scene.camera.view_up = [0.0, 0.0, 1.0]
    # scene.scene.camera.clipping_range = [1.5738238620907348, 6.861972426889951]
    # scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    mlab.savefig(
        "/Users/joonas/Documents/Manuscripts/ThermalNoise/figures/validation/cylinder.png",
        size=(800, 800),
    )

    Np = 30

    x = np.linspace(-0.95, 0.95, Np)
    fp = np.array((x, np.zeros(x.shape), np.zeros(x.shape))).T

    B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)
    B = noise_var(mesh, B_coupling, vl)

    # B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

    a = 0.5
    L = 2
    rat = L / (2 * a)
    Gfact = (
        1
        / (8 * np.pi)
        * (
            (3 * rat ** 5 + 5 * rat ** 3 + 2) / (rat ** 2 * (1 + rat ** 2) ** 2)
            + 3 * np.arctan(rat)
        )
    )
    Ban = np.sqrt(Gfact) * mu0 * np.sqrt(kB * T * sigma * d) / a

    plt.figure(figsize=(5, 5))
    plt.plot(x, Ban * np.ones(x.shape) * 1e15, label="Analytic", linewidth=2)
    plt.plot(
        x,
        np.sqrt(B[:, 0]) * 1e15,
        "x",
        label="Numerical",
        markersize=10,
        markeredgewidth=2,
    )
    plt.grid()
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    plt.legend(frameon=False)
    plt.xlabel("Distance along long axis")
    plt.ylabel("DC noise along axis (fT/rHz)")
    plt.tight_layout()

    plt.figure()
    plt.semilogy(x, np.sqrt(B[:, 0]), label="x")
    plt.semilogy(x, np.sqrt(B[:, 1]), label="y")
    plt.semilogy(x, np.sqrt(B[:, 2]), "--", label="z")
    plt.legend()
    plt.xlabel("Distance along long axis x")
    plt.ylabel("DC noise (T/rHz)")



Unit disc, AC mode
------------------


.. code-block:: default


    mesh = trimesh.load(
        pkg_resources.resource_filename(
            "bfieldtools", "example_meshes/unitdisc_extremelyfine.stl"
        )
    )


    # Nfreqs = 100
    # freqs = np.logspace(0, 3, Nfreqs) #30 frequencies from 1 to 1000 Hz
    # inds = np.where(freqs < 600)
    # freqs = freqs[inds]
    # Nfreqs = freqs.shape[0]

    Nfreqs = 70
    freqs = np.linspace(0, 1200, Nfreqs)

    S = np.ones(mesh.triangles_center.shape[0]) * sigma
    sheet_resistance = 1 / (d * S)

    # Compute the resistance and inductance matrices
    R = resistance_matrix(mesh, sheet_resistance=sheet_resistance)
    M = self_inductance_matrix(mesh, Nchunks=Nchunks, quad_degree=quad_degree)

    vl = compute_current_modes_ind_res(mesh, M, R, freqs, T, closed=False)

    #
    # fp = np.zeros((1,3))
    # fp[0,2] = 0.1

    Np = 20
    z = np.linspace(0.05, 0.2, Np)
    fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

    B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

    Bf = np.sqrt(noise_var(mesh, B_coupling, vl))

    # r = 1
    # Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*fp[0,2]**2))*(1/(1+fp[0,2]**2/r**2))

    plt.figure(figsize=(5, 5))
    plt.loglog(freqs, Bf[:, 2, :].T * 1e15, linewidth=2)
    plt.grid()
    plt.ylim(1, 20)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    plt.legend(frameon=False)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel(r"$B_z$ noise (fT/rHz)")
    plt.tight_layout()

    cutf = np.zeros(Np)
    for i in range(Np):
        idx = np.max(np.where(Bf[i, 2, :] >= 1 / np.sqrt(2) * Bf[i, 2, 0]))
        cutf[i] = freqs[idx]

    cutf_an = 1 / (4 * mu0 * sigma * d * z)

    plt.figure(figsize=(5, 5))
    plt.loglog(z, cutf_an, linewidth=2, label="Infinite plane")
    plt.loglog(z, cutf, "x", markersize=10, markeredgewidth=2, label="Disc")
    plt.grid()
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    plt.legend(frameon=False)
    plt.xlabel("Distance (z/R)")
    plt.ylabel("3-dB cutoff frequency (Hz)")
    plt.tight_layout()


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 19 minutes  14.276 seconds)


.. _sphx_glr_download_auto_examples_thermal_noise_thermal_noise_simulation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: thermal_noise_simulation.py <thermal_noise_simulation.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: thermal_noise_simulation.ipynb <thermal_noise_simulation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
