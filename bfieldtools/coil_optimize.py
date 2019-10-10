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


def optimize_streamfunctions(meshobj,
                             bfield_specification,
                             objective='minimum_inductive_energy',
                             solver=None,
                             solver_opts={},
                             problem=None,
                             boundary_constraints='all_zero'):
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
    solver: string
        string specifying which solver CVXPy will use
    solver_opt: dict
        dict containing solver options CVXPY will pass to the solver
    problem: CVXPY problem object
        If passed, will use already existing problem (MUST BE SAME DIMENSIONS) to
        skip DCP processing/reformulation time.
    boundary_constraints: string or dict
        Specifies how boundaries are handled. If 'all_zero' (default), boundary vertices are not included
        in the optimization, and are set at zero. If dict, must have entries
            'zero_eq_indices' (list of ints)
            'iso_eq_indices'  (list of lists of ints)
        which contain the vertex indices for which the stream function is set to zero or
        to be equal across elements.

    Returns
    -------
    I: vector
        Vector with length len(`meshobj.mesh.vertices`), containing the optimized current density values
        at each mesh vertex
    prob: CVXPY problem object
        CVXPY problem object containing data, formulation, solution, metric etc

    '''

    if objective == 'minimum_inductive_energy':
        objective = (1, 0)
    elif objective == 'minimum_resistive_energy':
        objective = (0, 1)

    #Initialize inequality constraint matrix and constraints
    constraint_matrix = np.zeros((0, len(meshobj.mesh.vertices)))
    upper_bounds = np.zeros((0, ))
    lower_bounds = np.zeros((0, ))

    #If boundaries are all zero, don't include them in the optimization
    if boundary_constraints == 'all_zero':
        indices = meshobj.inner_verts
    else:
        indices = np.arange(0, len(meshobj.mesh.vertices))

    #Populate inequality constraints with bfield specifications
    for spec in bfield_specification:


        C = spec['C'][:, indices]

        #Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
        #If not 3D matrix, assuming the use of spherical harmonics
        if C.ndim == 3:

            C = C.transpose((1, 0, 2))
            C = C.reshape((C.shape[0], -1)).T

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
        constraint_matrix = np.append(constraint_matrix, C, axis=0)

        upper_bounds = np.append(upper_bounds, upper_bound, axis=0)
        lower_bounds = np.append(lower_bounds, lower_bound, axis=0)


    #Construct quadratic objective matrix
    if objective == (1, 0):

        quadratic_matrix = meshobj.inductance[indices][:, indices]

    elif objective == (0, 1):

        quadratic_matrix = meshobj.resistance[indices][:, indices]

    else:

        L = meshobj.inductance[indices][:, indices]

        R = meshobj.resistance[indices][:, indices]

        print('Scaling inductance and resistance matrices before optimization. This requires eigenvalue computation, hold on.')

        max_eval_L = largest_eigh(L, eigvals=(L.shape[0]-1, L.shape[0]-1))[0][0]
        max_eval_R = largest_eigh(R, eigvals=(L.shape[0]-1, L.shape[0]-1))[0][0]

        scaled_R = max_eval_L / max_eval_R * R

        quadratic_matrix = (objective[0] * L  + objective[1] * scaled_R)


    #Scale whole quadratic term according to largest eigenvalue
    max_eval_quad = largest_eigh(quadratic_matrix, eigvals=(quadratic_matrix.shape[0]-1, quadratic_matrix.shape[0]-1))[0][0]

    quadratic_matrix /= max_eval_quad


    #Compute, scale constraint matrix according to largest singular value
    u, s, vt = svds(constraint_matrix, k=1)

    #If no pre-constructed problem is passed, create it
    if problem is None:
        print('Pre-existing problem not passed, creating...')
        #Symbolic variable for CVXPY
        x = cp.Variable(shape=(len(quadratic_matrix), ), name='x')

        #Parameters into which data is loaded
        P = cp.Parameter(shape=quadratic_matrix.shape, name='P', PSD=True)
        G = cp.Parameter(shape=constraint_matrix.shape, name='G')
        lb = cp.Parameter(shape=lower_bounds.shape, name ='lb')
        ub = cp.Parameter(shape=upper_bounds.shape, name='ub')

        #Formulate problem and constraints
        objective = cp.Minimize((1/2)*cp.quad_form(x, P))

        constraints = [G@x >= lb, G@x <= ub]

        if type(boundary_constraints) == dict:
            if boundary_constraints['zero_eq_indices'] is not None:
                print('Passing boundary zero equality constraint')
                constraints += [x[boundary_constraints['zero_eq_indices']] == np.zeros_like(boundary_constraints['zero_eq_indices'])]

            if boundary_constraints['iso_eq_indices'] is not None:
                print('Passing boundary equality constraints')
                for i in range(len(boundary_constraints['iso_eq_indices'])):
                    for j in range(1, len(boundary_constraints['iso_eq_indices'][i])):
                        constraints += [x[boundary_constraints['iso_eq_indices'][i][j]] == x[boundary_constraints['iso_eq_indices'][i][0]]]

        problem = cp.Problem(objective,
                          constraints)
    else:
        print('Existing problem passed')

    #Assign values to parameters
    print('Passing parameters to problem...')
    for par in problem.parameters():
        if par.name() == 'P':
            #Make sure that quadratic matrix is positive semi-definite, scale constraint matrix
            par.value = .5 * (quadratic_matrix + quadratic_matrix.T)
        elif par.name() == 'G':
            par.value = constraint_matrix / s[0]
        elif par.name() == 'lb':
            par.value = lower_bounds
        elif par.name() == 'ub':
            par.value = upper_bounds
        else:
            print('Unknown parameter encountered')


    #Run solver
    print('Passing problem to solver...')
    problem.solve(solver=solver, verbose=True, **solver_opts)

    #Build final I vector with zeros on boundary elements, scale by same singular value as constraint matrix

    I = problem.variables()[0].value / s[0]


    return I, problem