
import pytest
from bfieldtools import coil_optimize
from bfieldtools.tests.test_conductor import _fake_conductor
import cvxpy
import numpy as np

from numpy.testing import assert_allclose

def test_coil_optimize():
    """
    test that optimization works for different conductor objects and solvers
    """
    
    
    
    for test_mesh in ['unit_sphere', 'unit_disc', 'plane_w_holes']:
        results = []
        
        for basis_name in ['suh']:#, 'inner', 'vertex']:
            
            c = _fake_conductor(mesh_name=test_mesh, basis_name=basis_name, N_suh=10)
            
            spec = dict(coupling=c.B_coupling(np.array([[0, 0, 2]])), target=np.array([[0, 0, 1]]), abs_error=0.01)
            
            for objective in [(0, 1)]:#, (1, 0), (0.5, 0.5)]:
                
                #For now, test with all solvers that can handle SOC problems
                for solver in [i for i in cvxpy.solvers.defines.INSTALLED_CONIC_SOLVERS if i != 'GLPK' and i!= 'GLPK_MI']:
                    results.append(coil_optimize.optimize_streamfunctions(c,
                                                                        [spec],
                                                                        objective,
                                                                        solver
                                                                        )[0])
        if len(results) > 1:
            #tolerance is quite high, since some solvers give a bit differing results
            #in real life, let's not use those solvers.
            assert_allclose(results[-2], results[-1], rtol=2e-1)
    
