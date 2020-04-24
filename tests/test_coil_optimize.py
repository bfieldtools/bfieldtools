import pytest
from bfieldtools import coil_optimize
from .test_mesh_conductor import _fake_mesh_conductor
import cvxpy
import numpy as np

from numpy.testing import assert_allclose


def test_coil_optimize():
    """
    test that optimization works for different mesh_conductor objects and solvers
    """

    for test_mesh in ["unit_sphere", "unit_disc", "plane_w_holes"]:

        for basis_name in ["suh"]:  # , 'inner', 'vertex']:

            c = _fake_mesh_conductor(
                mesh_name=test_mesh, basis_name=basis_name, N_suh=10
            )

            spec = dict(
                coupling=c.B_coupling(np.array([[0, -0.1, 1], [0, 0, 1], [0, 0.1, 1]])),
                target=np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1]]),
                abs_error=0.01,
            )

            for objective in [
                "minimum_inductive_energy",
                "minimum_ohmic_power",
                (0.5, 0.5),
            ]:
                results = []

                results.append(coil_optimize.optimize_lsq(c, [dict(spec)], reg=1e3))

                # For now, test with all solvers that can handle SOC problems
                for solver in [
                    i
                    for i in cvxpy.solvers.defines.INSTALLED_CONIC_SOLVERS
                    if i not in ("GLPK", "GLPK_MI", "SCS")
                ]:
                    results.append(
                        coil_optimize.optimize_streamfunctions(
                            c, [spec], objective, solver
                        )[0]
                    )
                if len(results) > 1:
                    # tolerance is quite high, since some solvers give a bit differing results
                    # in real life, let's not use those solvers.
                    assert_allclose(
                        results[-2],
                        results[-1],
                        rtol=5e-1,
                        atol=0.001 * np.mean(np.abs(results[-2])),
                    )


def test_standalone_functions():
    """
    Tests standalone functions in coil_optimize

    """

    P = 2 * np.array([[2, 0.5], [0.5, 1]])
    q = np.array([1.0, 1.0])
    G = np.array([[-1.0, 0.0], [0.0, -1.0]])
    h = np.array([0.0, 0.0])
    A = np.array([[1.0, 1.0], [1, 2]])
    b = np.array([1.0, 0])

    coil_optimize.cvxopt_solve_qp(P, q)
    coil_optimize.cvxopt_solve_qp(P, q, G, h)
    coil_optimize.cvxpy_solve_qp(P, G, h, solver="CVXOPT")
