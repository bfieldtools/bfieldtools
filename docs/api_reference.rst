API Reference
=============

   


:mod:`bfieldtools` -- Main package
----------------------------------

.. contents::
   :local:
   :depth: 2

:py:mod:`bfieldtools`:

.. automodule:: bfieldtools
   :no-members:
   :no-inherited-members:

Surface mesh field calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.mesh_conductor`:



.. automodule:: bfieldtools.mesh_conductor
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.mesh_conductor

.. autosummary::
   :toctree: generated/

   MeshConductor
   CouplingMatrix
   StreamFunction
   
 
Convex optimization of stream functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.coil_optimize`:

.. automodule:: bfieldtools.coil_optimize
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.coil_optimize

.. autosummary::
   :toctree: generated/
   
   optimize_streamfunctions
   optimize_lsq
   cvxpy_solve_qp
   cvxopt_solve_qp
   


Mesh impedance
~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.mesh_impedance`:

.. automodule:: bfieldtools.mesh_impedance
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.mesh_impedance


.. autosummary::
   :toctree: generated/
   
   resistance_matrix
   self_inductance_matrix
   mutual_inductance_matrix
   mesh2line_mutual_inductance

   
   
Surface mesh calculus
~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.mesh_calculus`:



.. automodule:: bfieldtools.mesh_calculus
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.mesh_calculus

.. autosummary::
   :toctree: generated/

   laplacian_matrix
   mass_matrix
   gradient_matrix
   gradient
   divergence_matrix
   divergence
   curl_matrix
   curl
   


Stream function magnetic couplings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.mesh_magnetics`:

.. automodule:: bfieldtools.mesh_magnetics
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.mesh_magnetics

.. autosummary::
   :toctree: generated/
   
	magnetic_field_coupling
	magnetic_field_coupling_analytic
	scalar_potential_coupling
	vector_potential_coupling


Analytical integrals for fields and potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.integrals`:

.. automodule:: bfieldtools.integrals
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.integrals

.. autosummary::
   :toctree: generated/
   
	gamma0
	omega
	triangle_potential_uniform
	triangle_potential_approx
	potential_dipoles
	potential_vertex_dipoles
	triangle_potential_dipole_linear


Scalar function contouring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.contour`:

.. automodule:: bfieldtools.contour
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.contour

.. autosummary::
   :toctree: generated/
   
   scalar_contour
   simplify_contour
   
   
Line currents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.line_conductor`:

.. automodule:: bfieldtools.line_conductor
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/
	
	LineConductor

:py:mod:`bfieldtools.line_magnetics`:

.. automodule:: bfieldtools.line_magnetics
   :no-members:
   :no-inherited-members:

.. currentmodule:: bfieldtools.line_magnetics

.. autosummary::
   :toctree: generated/
   
Spherical harmonics
~~~~~~~~~~~~~~~~~~~


:py:mod:`bfieldtools.sphtools`:

.. automodule:: bfieldtools.sphtools
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/
	
	SphBasis
	cartesian2spherical
	spherical2cartesian
	Rotmatrix
	cartvec2sph
	sphvec2cart
	lpmn_em
	derlpmn_em
	xlm
	ylm
	derxlm
	sinxlm
	dthylm
	dphiylm
	Plm
	Blm
	Wlm
	Vlm
	potential
	field
	compute_sphcoeffs_mesh
	basis_potentials
	basis_fields
	plotYlms
	plotYlm
	plotWlm
	plotBWlm_volume
	plotVlm
	plotBVlm_volume
	
	
Surface harmonics
~~~~~~~~~~~~~~~~~


:py:mod:`bfieldtools.suhtools`:

.. automodule:: bfieldtools.suhtools
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/
   
   SuhBasis
   
   
Thermal noise
~~~~~~~~~~~~~~~~~


:py:mod:`bfieldtools.thermal_noise`:

.. automodule:: bfieldtools.thermal_noise
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/
   
   compute_current_modes
   noise_covar
   noise_var
   noise_covar_dir
   visualize_current_modes
   
   
Visualization and plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.viz`:



.. automodule:: bfieldtools.viz
   :no-members:
   :no-inherited-members:
   
.. currentmodule:: bfieldtools.viz
   
.. autosummary::
   :toctree: generated/
   
   plot_mesh
   plot_3d_current_loops
   plot_cross_section
   plot_data_on_vertices
   plot_data_on_faces
   
Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:mod:`bfieldtools.utils`:



.. automodule:: bfieldtools.utils
   :no-members:
   :no-inherited-members:
   
.. currentmodule:: bfieldtools.utils
   
.. autosummary::
   :toctree: generated/
   
   combine_meshes
   get_quad_points
   get_line_quad_points
   dual_areas
   find_mesh_boundaries
   load_example_mesh
   fix_normals
   inner2vert
   vert2inner
   MeshProjection
