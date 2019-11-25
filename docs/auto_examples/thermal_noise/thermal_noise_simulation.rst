.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_thermal_noise_thermal_noise_simulation.py>` to download the full example code
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

    from bfieldtools.thermal_noise import compute_current_modes, visualize_current_modes, compute_dc_Bnoise, compute_ac_Bnoise

    import pkg_resources



    #Fix the simulation parameters
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4*np.pi*1e-7








Unit sphere
------------


.. code-block:: default



    Np = 10
    R = np.linspace(0.1, 1, Np)
    fp = np.zeros((1,3))

    B = np.zeros((Np,3))
    for i in range(Np):
        mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_sphere.stl'))
        mesh.apply_scale(R[i])
        vl = compute_current_modes(mesh)

        vl[:,0] = np.zeros(vl[:,0].shape) # fix DC-component


        Btemp = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)
        B[i] = Btemp

    visualize_current_modes(mesh,vl, 40, 5)

    Ban = mu0*np.sqrt(2*sigma*d*kB*T/(3*np.pi*(R)**2))

    plt.figure()
    plt.semilogy(R, Ban,label='Analytic')
    plt.semilogy(R, B[:,2],'x',label='Numerical')
    plt.legend()
    plt.xlabel('Sphere radius')
    plt.ylabel('DC noise Bz (T/rHz)')
    plt.tight_layout()


    RE = np.abs((B[:,2]-Ban))/np.abs(Ban)*100
    plt.figure()
    plt.plot(R, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)
    plt.xlabel('Sphere radius')
    plt.ylabel('Relative error (%)')




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

    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.



Unit disc, DC noise
---------------------


.. code-block:: default


    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_disc.stl'))
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    vl = compute_current_modes(mesh)

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    visualize_current_modes(mesh,vl, 40, 5)

    Np = 30

    z = np.linspace(0.1, 1, Np)
    fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

    B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

    r = 1
    Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*z**2))*(1/(1+z**2/r**2))

    plt.figure()
    plt.semilogy(z, Ban,label='Analytic')
    plt.semilogy(z, B[:,2],'x',label='Numerical')
    plt.legend()
    plt.xlabel('Distance d/R')
    plt.ylabel('DC noise Bz (T/rHz)')
    plt.tight_layout()

    plt.figure()
    plt.plot(z, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)
    plt.xlabel('Distance d/R')
    plt.ylabel('Relative error (%)')




.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_004.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_005.png
            :class: sphx-glr-multi-img

.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_006.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1207 vertices by 30 target points... took 0.04 seconds.



Closed cylinder, DC noise
---------------------


.. code-block:: default


    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder.stl'))
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    vl = compute_current_modes(mesh)

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    visualize_current_modes(mesh,vl, 8, 1)

    Np = 30

    x = np.linspace(-0.95, 0.95, Np)
    fp = np.array((x,np.zeros(x.shape), np.zeros(x.shape))).T

    B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

    a = 0.5
    L = 2
    rat = L/(2*a)
    Gfact = 1/(8*np.pi) * ((3*rat**5+5*rat**3+2)/(rat**2*(1+rat**2)**2) + 3*np.arctan(rat))
    Ban = np.sqrt(Gfact)*mu0*np.sqrt(kB*T*sigma*d)/a

    plt.figure()
    plt.semilogy(x, Ban*np.ones(x.shape),label='Analytic',linewidth = 2)
    plt.semilogy(x, B[:,0],'x',label='Numerical')
    plt.legend()
    plt.xlabel('Distance along long axis')
    plt.ylabel('DC noise long axis (T/rHz)')
    plt.tight_layout()

    plt.figure()
    plt.semilogy(x, B[:,0],label='x')
    plt.semilogy(x, B[:,1],label='y')
    plt.semilogy(x, B[:,2],'--',label='z')
    plt.legend()
    plt.xlabel('Distance along long axis x')
    plt.ylabel('DC noise (T/rHz)')






.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_007.png
            :class: sphx-glr-multi-img

    *

      .. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_008.png
            :class: sphx-glr-multi-img

.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_009.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    face_normals didn't match triangles, ignoring!
    Computing magnetic field coupling matrix, 3842 vertices by 30 target points... took 0.16 seconds.



Unit disc, AC mode
------------------


.. code-block:: default


    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_disc.stl'))
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)


    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))
    vl = compute_current_modes(mesh)

    fp = np.zeros((1,3))
    fp[0,2] = 0.1

    Nfreqs = 30
    freqs = np.logspace(0, 3, Nfreqs) #30 frequencies from 1 to 1000 Hz

    Bf = compute_ac_Bnoise(mesh,vl,fp,freqs,sigma,d,T)

    r = 1
    Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*fp[0,2]**2))*(1/(1+fp[0,2]**2/r**2))

    plt.figure()
    plt.loglog(freqs,Bf[:,0,2],label = 'Numerical')
    plt.loglog(freqs, Ban*np.ones(freqs.shape), '--',label = 'Analytical, DC')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Field noise (T/rHz)')
    plt.legend()
    plt.grid(which='both')
    plt.tight_layout()


.. code-block:: pytb

    Traceback (most recent call last):
      File "/l/conda-envs/mne/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 474, in _memory_usage
        multiprocess=True)
      File "/l/conda-envs/mne/lib/python3.6/site-packages/memory_profiler.py", line 336, in memory_usage
        returned = f(*args, **kw)
      File "/l/conda-envs/mne/lib/python3.6/site-packages/sphinx_gallery/gen_rst.py", line 465, in __call__
        exec(self.code, self.globals)
      File "/l/bfieldtools/examples/thermal_noise/thermal_noise_simulation.py", line 174, in <module>
        Bf = compute_ac_Bnoise(mesh,vl,fp,freqs,sigma,d,T)
      File "/l/bfieldtools/bfieldtools/thermal_noise.py", line 255, in compute_ac_Bnoise
        B[j, :, 0] += (B_coupling[:, :, 0] @vec)**2
    ValueError: shapes (1,3) and (1207,) not aligned: 3 (dim 1) != 1207 (dim 0)





.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  20.474 seconds)

**Estimated memory usage:**  692 MB


.. _sphx_glr_download_auto_examples_thermal_noise_thermal_noise_simulation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: thermal_noise_simulation.py <thermal_noise_simulation.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: thermal_noise_simulation.ipynb <thermal_noise_simulation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
