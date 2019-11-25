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






.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    SVG path loading unavailable!
    Traceback (most recent call last):
      File "/u/76/zetterr1/unix/.local/lib/python3.6/site-packages/trimesh/path/exchange/svg_io.py", line 18, in <module>
        from svg.path import parse_path
    ModuleNotFoundError: No module named 'svg'



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

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

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
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.07 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.07 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.07 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.07 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.07 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    /l/bfieldtools/bfieldtools/thermal_noise.py:69: RuntimeWarning: invalid value encountered in sqrt
      vl[inner_verts, i] = v[:, i]/np.sqrt(u[i])
    Computing magnetic field coupling matrix, 2562 vertices by 1 target points... took 0.08 seconds.
    0 0
    1 0
    2 0
    3 0
    4 0
    5 0
    6 0
    0 1
    1 1
    2 1
    3 1
    4 1
    5 1
    6 1
    0 2
    1 2
    2 2
    3 2
    4 2
    5 2
    6 2
    0 3
    1 3
    2 3
    3 3
    4 3
    5 3
    6 3
    0 4
    1 4
    2 4
    3 4
    4 4
    5 4
    6 4
    0 5
    1 5
    2 5
    3 5
    4 5



Unit disc, DC noise
---------------------


.. code-block:: default


    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_disc.stl'))
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    vl = compute_current_modes(mesh)

    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                   size=(800, 800))

    visualize_current_modes(mesh,vl, 42, 5, contours=False)

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

    0 0
    1 0
    2 0
    3 0
    4 0
    5 0
    6 0
    0 1
    1 1
    2 1
    3 1
    4 1
    5 1
    6 1
    0 2
    1 2
    2 2
    3 2
    4 2
    5 2
    6 2
    0 3
    1 3
    2 3
    3 3
    4 3
    5 3
    6 3
    0 4
    1 4
    2 4
    3 4
    4 4
    5 4
    6 4
    0 5
    1 5
    2 5
    3 5
    4 5
    5 5
    6 5
    Computing magnetic field coupling matrix, 1207 vertices by 30 target points... took 0.04 seconds.



Closed cylinder, DC noise
--------------------------


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
    0 0
    1 0
    2 0
    0 1
    1 1
    2 1
    0 2
    1 2
    Computing magnetic field coupling matrix, 3842 vertices by 30 target points... took 0.14 seconds.



Unit disc, AC mode
------------------


.. code-block:: default


    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_disc.stl'))
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)


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


.. image:: /auto_examples/thermal_noise/images/sphx_glr_thermal_noise_simulation_010.png
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Computing magnetic field coupling matrix, 1207 vertices by 1 target points... took 0.06 seconds.
    Calculating potentials, chunk 1/1
    Frequency 1.000000 computed
    Frequency 1.268961 computed
    Frequency 1.610262 computed
    Frequency 2.043360 computed
    Frequency 2.592944 computed
    Frequency 3.290345 computed
    Frequency 4.175319 computed
    Frequency 5.298317 computed
    Frequency 6.723358 computed
    Frequency 8.531679 computed
    Frequency 10.826367 computed
    Frequency 13.738238 computed
    Frequency 17.433288 computed
    Frequency 22.122163 computed
    Frequency 28.072162 computed
    Frequency 35.622479 computed
    Frequency 45.203537 computed
    Frequency 57.361525 computed
    Frequency 72.789538 computed
    Frequency 92.367086 computed
    Frequency 117.210230 computed
    Frequency 148.735211 computed
    Frequency 188.739182 computed
    Frequency 239.502662 computed
    Frequency 303.919538 computed
    Frequency 385.662042 computed
    Frequency 489.390092 computed
    Frequency 621.016942 computed
    Frequency 788.046282 computed
    Frequency 1000.000000 computed




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 1 minutes  32.243 seconds)

**Estimated memory usage:**  2312 MB


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
