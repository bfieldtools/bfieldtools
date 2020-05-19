"""
Includes files for coil optimization (stream function optimization)
using either a numerical solver or regularized least squares

"""

__all__ = [
    "cvxopt_solve_qp",
    "cvxpy_solve_qp",
    "optimize_lsq",
    "optimize_streamfunctions",
]

import numpy as np
import cvxopt
from cvxopt import matrix
from scipy.sparse.linalg import svds
from scipy.linalg import eigh
from scipy.linalg import eigvalsh
import cvxpy as cp

from .mesh_conductor import StreamFunction


def cvxpy_solve_qp(P, G, h, solver=cp.MOSEK, tolerance=None):
    """
    Bare-bones quadratic programming solver function for CVXPY.
    
    Minimizes

    .. math::
        (1/2) \mathbf{x}^T  \mathbf{P}  \mathbf{x}

    subject to

    .. math::
        \mathbf{G} \mathbf{x} \leq \mathbf{h}


    Parameters
    ----------

    P: (N_x, N_x) array
        Quadratic minimization matrix
    G: (N_h, N_x) array
        Linear inequality matrix
    h: (N_h, ) array
        Linear inequality constraint vector
    solver: string or cvxpy solver
        solver to use in CVXPY, can be passed as string or directly
    tolerance: float or None (default)
        Override default tolerance values

    Returns
    -------
    x.value: (N_x, ) array
        Optimized solution for x (if available)
    prob: CVXPY problem
        Problem object with optimization info

    """

    P = 0.5 * (P + P.T)

    x = cp.Variable(len(P))

    objective = cp.Minimize((1 / 2) * cp.quad_form(x, P))

    constraints = [G @ x <= h]

    prob = cp.Problem(objective, constraints)
    if tolerance is None:
        prob.solve(solver=solver, verbose=True)  # Returns the optimal value.
    elif solver == cp.CVXOPT:
        prob.solve(
            solver=solver,
            verbose=True,
            abstol=tolerance,
            feastol=tolerance,
            reltol=tolerance,
        )  # Returns the optimal value.
    elif solver == cp.SCS:
        prob.solve(
            solver=solver, verbose=True, eps=tolerance
        )  # Returns the optimal value.

    # Print result.
    print("\nThe optimal value is", prob.value)
    print("A solution x is")
    print(x.value)
    print("A dual solution corresponding to the inequality constraints is")
    print(prob.constraints[0].dual_value)

    return x.value, prob


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, tolerance=1e-7, **kwargs):
    """
    Use cvxopt (without CVXPY wrapper) for quadratic programming.
    
    Minimize

    .. math::
        (1/2) \mathbf{x}^T  \mathbf{P}  \mathbf{x} + \mathbf{q}^T \mathbf{x}

    subject to

    .. math::
        \mathbf{G} \mathbf{x} \leq \mathbf{h}

    and

    .. math::
        \mathbf{A} \mathbf{x} = \mathbf{b}

    Parameters
    ----------

    P: (N_x, N_x) array
        Quadratic minimization matrix
    q: (N_x, ) array
        Linear penalty term vector
    G: (N_h, N_x) array
        Linear inequality matrix
    h: (N_h, ) array
        Linear inequality constraint vector
    A: (N_b, N_x) array
        Linear equality matrix
    b: (N_b, ) array
        Linear equality constraint vector
    **kwargs
        Use to pass cvxopt solver options, such as tolerance

    Returns
    -------
    x: (N_x, ) array
        Optimized solution for x (if available)
    sol: CVXOPT solution
        Solution object with optimization info

    """

    P = 0.5 * (P + P.T)  # make sure P is symmetric

    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])

    for key, val in kwargs:
        cvxopt.solver.options[key] = val

    sol = cvxopt.solvers.qp(*args)
    if "optimal" not in sol["status"]:
        return None, sol
    return np.array(sol["x"]).reshape((P.shape[1],)), sol


def optimize_streamfunctions(
    mesh_conductor,
    bfield_specification,
    objective="minimum_inductive_energy",
    solver=None,
    solver_opts={},
    problem=None,
):
    """
    Quadratic optimization of coil stream function according to a specified objective.

    Utilizes CVXPY and a numerical iterative solver.

    Parameters
    ----------
    mesh_conductor: MeshConductor object
        Contains Trimesh mesh as well as physical properties, e.g. inductance
    bfield_specification: list
        List in which element is a dictionary containing a coil specification.
        See notes for specification syntax.
    objective: string or dict
        if string, either *'minimum_inductive_energy'* or *'minimum_ohmic_power'*
        if tuple, should contain: (a, b), where a and b are floats describing the
        inductive and resitive weighting factors.
        The resistance matrix is scaled according to the largest singular value
        of the inductance matrix for consistent behavior across meshes.
    solver: string
        string specifying which solver CVXPY will use
    solver_opt: dict
        dict containing solver options CVXPY will pass to the solver
    problem: CVXPY problem object
        If passed, will use already existing problem (**MUST BE SAME DIMENSIONS**) to
        skip DCP processing/reformulation time.

    Returns
    -------
    s: vector
        Vector with length len(`mesh_conductor.mesh.vertices`), containing the
        optimized current density values at each mesh vertex
    prob: CVXPY problem object
        CVXPY problem object containing data, formulation, solution, metric etc

    Notes
    -----

    Each specification is a dict, which contains
     - coupling: Coupling matrix (N_r, N_verts, 3)
     - target: (N_r, 3)
     - abs_error: float or (N_r, 3)
     - rel_error: float or (N_r, 3)
    Either abs_error, rel_error or both must be present.

    """

    if objective == "minimum_inductive_energy":
        objective = (1, 0)
    elif objective == "minimum_ohmic_power":
        objective = (0, 1)

    quadratic_matrix = _construct_quadratic_objective(objective, mesh_conductor)

    constraint_matrix, upper_bounds, lower_bounds = _construct_constraints(
        mesh_conductor, bfield_specification
    )
    # Compute, scale constraint matrix according to largest singular value
    u, s, vt = svds(constraint_matrix, k=1)

    # If no pre-constructed problem is passed, create it
    if problem is None:
        print("Pre-existing problem not passed, creating...")
        # Symbolic variable for CVXPY
        x = cp.Variable(shape=(len(quadratic_matrix),), name="x")

        # Parameters into which data is loaded
        # P does not need to be PSD when in this formulation
        # It should be full rank, however
        P = cp.Parameter(shape=quadratic_matrix.shape, name="P")  # , PSD=True)
        G = cp.Parameter(shape=constraint_matrix.shape, name="G")
        lb = cp.Parameter(shape=lower_bounds.shape, name="lb")
        ub = cp.Parameter(shape=upper_bounds.shape, name="ub")

        # Formulate problem and constraints
        objective = cp.Minimize((1 / 2) * cp.sum_squares(P @ x))

        constraints = [G @ x >= lb, G @ x <= ub]

        problem = cp.Problem(objective, constraints)
    else:
        print("Existing problem passed")

    # Assign values to parameters
    print("Passing parameters to problem...")
    for par in problem.parameters():
        if par.name() == "P":
            # Make sure that quadratic matrix is positive semi-definite, scale constraint matrix
            # par.value = sqrtm(0.5 * (quadratic_matrix + quadratic_matrix.T))
            PP = 0.5 * (quadratic_matrix + quadratic_matrix.T)
            par.value = np.linalg.cholesky(PP).T
        elif par.name() == "G":
            par.value = constraint_matrix / s[0]
        elif par.name() == "lb":
            par.value = lower_bounds
        elif par.name() == "ub":
            par.value = upper_bounds
        else:
            print("Unknown parameter encountered")

    # Run solver
    print("Passing problem to solver...")
    problem.solve(solver=solver, verbose=True, **solver_opts)

    # extract optimized streamfunction, scale by same singular value as constraint matrix
    S = StreamFunction(
        problem.variables()[0].value / s[0], mesh_conductor=mesh_conductor
    )

    return S, problem


def optimize_lsq(
    mesh_conductor, bfield_specification, reg=1e3, objective="minimum_inductive_energy"
):
    """
    Optimization of coil stream function according to a specified objective using least-squares.
    
    Parameters
    ----------
    mesh_conductor: MeshConductor object
        Contains Trimesh mesh as well as physical properties, e.g. inductance
    bfield_specification: list
        List in which element is a dictionary containing a coil specification.
    objective: string or dict
        if string, either *'minimum_inductive_energy'* or *'minimum_ohmic_power'*
        if tuple, should contain: (a, b), where a and b are floats describing the
        inductive and resitive weighting factors.
        The resistance matrix is scaled according to the largest singular value
        of the inductance matrix for consistent behavior across meshes.
    reg: float
        Regularization/tradeoff parameter (lambda). A larger lambda leads to
        more emphasis on the specification, at the cost of the quadratic objective.
        The lambda value is relative to the maximum singular value w
        of the eigenvalue equation

        .. math::
            \mathbf{C}^T \mathbf{C}  \mathbf{v}[:,i] = \mathbf{w}[i]  \mathbf{Q} \mathbf{v}[:,i]

        where C is the constraint matrix and Q is the quadratic objective matrix.

    Returns
    -------
    S: StreamFunction
        Optimization solution

    Notes
    -----

    Each specification is a dict, which contains
     - coupling: Coupling matrix (N_r, N_verts, 3)
     - target: (N_r, 3)

    **NOTE** The following spec parameters are ignored:
     - abs_error: float or (N_r, 3)
     - rel_error: float or (N_r, 3)
    """

    if objective == "minimum_inductive_energy":
        objective = (1, 0)
    elif objective == "minimum_ohmic_power":
        objective = (0, 1)

    quadratic_matrix = _construct_quadratic_objective(objective, mesh_conductor)

    # Make sure that quadratic matrix is positive semi-definite
    quadratic_matrix = 0.5 * (quadratic_matrix + quadratic_matrix.T)

    print("Error tolerances in specification will be ignored when using lsq")
    for spec in bfield_specification:
        spec.pop("abs_error", None)
        spec.pop("rel_error", None)

    constraint_matrix, target = _construct_constraints(
        mesh_conductor, bfield_specification
    )

    # Compute, scale constraint matrix according to largest singular value
    u, s, vt = svds(constraint_matrix, k=1)
    constraint_matrix /= s[0]

    ss_max = eigvalsh(
        constraint_matrix.T @ constraint_matrix,
        quadratic_matrix,
        eigvals=[quadratic_matrix.shape[1] - 1, quadratic_matrix.shape[1] - 1],
    )

    S = np.linalg.solve(
        constraint_matrix.T @ constraint_matrix + ss_max / reg * quadratic_matrix,
        constraint_matrix.T @ target,
    )

    return StreamFunction(S / s[0], mesh_conductor)


def _construct_constraints(mesh_conductor, bfield_specification):
    """
    Stacks together constraint and coupling matrices for stream function optimization

    """
    # Initialize inequality constraint matrix and constraints
    constraint_matrix = np.zeros((0, mesh_conductor.basis.shape[1]))
    upper_bound_stack = np.zeros((0,))
    lower_bound_stack = np.zeros((0,))

    target_stack = np.zeros((0,))

    # Populate inequality constraints with bfield specifications
    for spec in bfield_specification:

        # Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
        # If not 3D matrix, assuming the use of spherical harmonics
        if spec["coupling"].ndim == 3:
            C = spec["coupling"].transpose((2, 0, 1))
            C = C.reshape((C.shape[0], -1)).T
        else:
            C = spec["coupling"]

        # Append specification to constraint matrix and bounds
        constraint_matrix = np.append(constraint_matrix, C, axis=0)

        if "rel_error" in spec or "abs_error" in spec:
            # Apply relative error to bounds
            if "rel_error" in spec:
                upper_bound = spec["target"] * (
                    1 + np.sign(spec["target"]) * spec["rel_error"]
                )
                lower_bound = spec["target"] * (
                    1 - np.sign(spec["target"]) * spec["rel_error"]
                )

                # If present, also apply absolute error
                if "abs_error" in spec:
                    upper_bound += spec["abs_error"]
                    lower_bound -= spec["abs_error"]

            # Apply absolute error to bounds
            else:
                upper_bound = spec["target"] + spec["abs_error"]
                lower_bound = spec["target"] - spec["abs_error"]

            # Flatten to match C matrix
            upper_bound = upper_bound.flatten()
            lower_bound = lower_bound.flatten()

            # Append to stack
            upper_bound_stack = np.append(upper_bound_stack, upper_bound, axis=0)
            lower_bound_stack = np.append(lower_bound_stack, lower_bound, axis=0)

        else:
            target = spec["target"]

            # Flatten to match C matrix
            target = target.flatten()

            # Append to stack
            target_stack = np.append(target_stack, target, axis=0)

    if "rel_error" in spec or "abs_error" in spec:
        return constraint_matrix, upper_bound_stack, lower_bound_stack

    return constraint_matrix, target_stack


def _construct_quadratic_objective(objective, mesh_conductor, deflate=True):
    """

    """
    # Construct quadratic objective matrix
    if objective == (1, 0):

        quadratic_matrix = mesh_conductor.inductance

    elif objective == (0, 1):

        quadratic_matrix = mesh_conductor.resistance

    elif isinstance(objective, tuple):

        L = mesh_conductor.inductance

        R = mesh_conductor.resistance

        print(
            "Scaling inductance and resistance matrices before optimization.\
              This requires eigenvalue computation, hold on."
        )

        max_eval_L = eigh(L, eigvals=(L.shape[0] - 1, L.shape[0] - 1))[0][0]
        max_eval_R = eigh(R, eigvals=(L.shape[0] - 1, L.shape[0] - 1))[0][0]

        scaled_R = max_eval_L / max_eval_R * R

        quadratic_matrix = objective[0] * L + objective[1] * scaled_R
    else:
        print("Custom objective passed, assuming it is a matrix of correct dimensions")
        quadratic_matrix = objective

    # Scale whole quadratic term according to largest eigenvalue
    max_eval_quad = eigh(
        quadratic_matrix,
        eigvals=(quadratic_matrix.shape[0] - 1, quadratic_matrix.shape[0] - 1),
    )[0][0]

    quadratic_matrix /= max_eval_quad

    if deflate:
        evals = eigh(quadratic_matrix, eigvals=(0, 1), eigvals_only=True)
        if abs(evals[0]) < 1e-8 * abs(evals[1]):
            quadratic_matrix += np.ones / np.sqrt(quadratic_matrix.shape[0])

    return quadratic_matrix
