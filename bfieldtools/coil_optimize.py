import numpy as np
import cvxopt
from cvxopt import matrix
from scipy.sparse.linalg import svds
from scipy.linalg import eigh as largest_eigh
import quadprog
import cvxpy as cp

def cvxpy_solve_qp(P, G, h, solver=cp.MOSEK, tolerance=None):

    P = .5 * (P + P.T)

    x = cp.Variable(len(P))

    objective = cp.Minimize((1/2)*cp.quad_form(x, P))

    constraints = [G@x <= h]

    prob = cp.Problem(objective,
                      constraints)
    if tolerance is None:
        prob.solve(solver=solver, verbose=True)  # Returns the optimal value.
    elif solver==cp.CVXOPT:
        prob.solve(solver=solver, verbose=True, abstol=tolerance, feastol=tolerance, reltol=tolerance)  # Returns the optimal value.
    elif solver==cp.SCS:
        prob.solve(solver=solver, verbose=True, eps=tolerance)  # Returns the optimal value.


    # Print result.
    print("\nThe optimal value is", prob.value)
    print("A solution x is")
    print(x.value)
    print("A dual solution corresponding to the inequality constraints is")
    print(prob.constraints[0].dual_value)

    return x.value, prob


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, sw=None, reg=None, tolerance=1e-7):
    '''
    Use cvxopt to minimize
    (1/2) * x' * P * x + q' * x

    subject to
    G * x <= h

    and
    A * x = b
    '''

    P = .5 * (P + P.T)  # make sure P is symmetric

    if sw is not None:
        n, m = P.shape[0], G.shape[0]

        E, Z = np.eye(m), np.zeros((m, n))

        P = np.vstack([np.hstack([P, Z.T]), np.hstack([Z, reg * np.eye(m)])])
        q = np.hstack([q, -sw * np.ones(m)])

        G = np.hstack([Z, E])
        h = np.zeros(m)

        A = np.hstack([G, -E])
        b = h


    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])

#    #For now, use rough tolerance setting, i.e. just set all arguments to be the same
    cvxopt.solvers.options['abstol'] = tolerance
    cvxopt.solvers.options['feastol'] = tolerance * 1e3
    cvxopt.solvers.options['reltol'] = tolerance

    #cvxopt.solvers.options['maxiters'] = 1000
    #cvxopt.solvers.options['kktreg'] = 1e-8


    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None, sol
    return np.array(sol['x']).reshape((P.shape[1],)), sol


def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    '''
    Use quadprog to minimize
    (1/2) * x' * P * x + q' * x

    subject to
    G * x <= h

    and
    A * x = b
    '''
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        if A.ndim == 1:
            A = A.reshape((1, A.shape[0]))
        if G is None:
            qp_C = -A.T
            qp_b = -b
        else:
            qp_C = -np.vstack([A, G]).T
            qp_b = -np.hstack([b, h])
        meq = A.shape[0]

    else:  # no equality constraint
        qp_C = -G.T if G is not None else None
        qp_b = -h if h is not None else None
        meq = 0

    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]


def optimize_streamfunctions(meshobj, bfield_specification,
                             objective='minimum_inductive_energy',
                             solver=None,
                             solver_opts={}):
    '''
    Quadratic optimization of coil stream function according to a specified objective,
    while keeping specified target field at target points within given constraints.

    Parameters
    ----------
    meshobj: MeshWrapper object
        Contains Trimesh mesh
    bfield_specification: list
        List in which element is a dictionary containing a field specification.
        Each dict contains:
        C: Coupling matrix (N_r, N_verts, 3)
        target_field: (N_r, 3)
        abs_error: float or (N_r, 3)
        rel_error: float or (N_r, 3)
    objective: string or dict
        if string, either 'minimum_inductive_energy' or 'minimum_resistive_energy'
        if tuple, should contain: (a, b), where a and b are floats describing the inductive and resitive weighting factors.
        The resistance matrix is scaled according to the largest singular value of the inductance matrix for consistent behavior
        across meshes.
    solver
    tolerance: float

    Returns
    -------
    I: vector
        Vector with length len(`meshobj.mesh.vertices`), containing the optimized current density values
        at each mesh vertex
    sol: dict
        Dict containing solution info and diagnostics supplied by cvxopt

    '''

    if objective == 'minimum_inductive_energy':
        objective = (1, 0)
    elif objective == 'minimum_resistive_energy':
        objective = (0, 1)

    #Initialize inequality constraint matrix and constraints
    constraint_matrix = np.zeros((0, len(meshobj.inner_verts)))
    upper_bounds = np.zeros((0, ))
    lower_bounds = np.zeros((0, ))


    #Populate inequality constraints with bfield specifications
    for spec in bfield_specification:

        #Limit C matrix to inner vertices (boundaries are kept at zero)
        inner_C = spec['C'][:, meshobj.inner_verts]

        #Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
        #If not 3D matrix, assuming the use of spherical harmonics
        if inner_C.ndim == 3:

            inner_C = inner_C.transpose((1, 0, 2))
            inner_C = inner_C.reshape((inner_C.shape[0], -1)).T

        #Apply relative error to bounds
        if spec['rel_error'] is not None:
            upper_bound = spec['target_field'] * (1 + np.sign(spec['target_field']) * spec['rel_error'])
            lower_bound = spec['target_field'] * (1 - np.sign(spec['target_field']) * spec['rel_error'])

        #Apply absolute error to bounds
        if spec['abs_error'] is not None:
            upper_bound = spec['target_field'] + spec['abs_error']
            lower_bound = spec['target_field'] - spec['abs_error']


        #Flatten to match C matrix
        upper_bound = upper_bound.flatten()
        lower_bound = lower_bound.flatten()


        # Append specification to constraint matrix and bounds
        constraint_matrix = np.append(constraint_matrix, inner_C, axis=0)

        upper_bounds = np.append(upper_bounds, upper_bound, axis=0)
        lower_bounds = np.append(lower_bounds, lower_bound, axis=0)


    #Construct quadratic objective matrix
    if objective == (1, 0):
        #Limit L matrix to inner vertices
        inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

        quadratic_matrix = inner_L

    elif objective == (0, 1):
        #Limit R matrix to inner vertices
        inner_R = meshobj.resistance[meshobj.inner_verts][:, meshobj.inner_verts]

        quadratic_matrix = inner_R

    else:
        #Limit L matrix to inner vertices
        inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

        #Limit R matrix to inner vertices
        inner_R = meshobj.resistance[meshobj.inner_verts][:, meshobj.inner_verts]

        print('Scaling inductance and resistance matrices before optimization. This requires eigenvalue computation, hold on.')

        max_eval_L = largest_eigh(inner_L, eigvals=(inner_L.shape[0]-1, inner_L.shape[0]-1))[0][0]
        max_eval_R = largest_eigh(inner_R, eigvals=(inner_L.shape[0]-1, inner_L.shape[0]-1))[0][0]

        scaled_R = max_eval_L / max_eval_R * inner_R

        quadratic_matrix = (objective[0] * inner_L  + objective[1] * scaled_R)


    #Scale whole quadratic term according to largest eigenvalue
    max_eval_quad = largest_eigh(quadratic_matrix, eigvals=(quadratic_matrix.shape[0]-1, quadratic_matrix.shape[0]-1))[0][0]

    quadratic_matrix /= max_eval_quad


    #Compute, scale constraint matrix according to largest singular value
    u, s, vt = svds(constraint_matrix, k=1)

    #Make sure that quadratic matrix is positive semi-definite, scale constraint matrix
    P = .5 * (quadratic_matrix + quadratic_matrix.T)
    G = constraint_matrix / s[0]

    #Symbolic variable for CVXPY
    x = cp.Variable(len(P))

    #Formulate problem and constraints
    objective = cp.Minimize((1/2)*cp.quad_form(x, P))


    constraints = [G@x >=lower_bounds, G@x <= upper_bounds]

    prob = cp.Problem(objective,
                      constraints)

#    prob.solve(solver=cp.OSQP, verbose=True, linsys_solver='mkl pardiso', rho=0.1, max_iter=50000)  # Returns the optimal value.

    #Run solver
    prob.solve(solver=solver, verbose=True, **solver_opts)

    #Build final I vector with zeros on boundary elements, scale by same singular value as constraint matrix
    I = np.zeros((meshobj.inductance.shape[0], ))
    I[meshobj.inner_verts] = x.value / s[0]


    return I, prob