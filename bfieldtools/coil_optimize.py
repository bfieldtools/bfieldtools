import numpy as np
import cvxopt
from cvxopt import matrix
from scipy.sparse.linalg import svds


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, tolerance=1e-7):
    '''
    Minimize
    (1/2) * x' * P * x + q' * x

    subject to
    G * x <= h

    and
    A * x = b
    '''

    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])

    #For now, use rough tolerance setting, i.e. just set all arguments to be the same
    cvxopt.solvers.options['abstol'] = tolerance
    cvxopt.solvers.options['feastol'] = tolerance
    cvxopt.solvers.options['reltol'] = tolerance

    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],)), sol

def optimize_streamfunctions(meshobj, bfield_specification,
                             objective='minimum_inductive_energy',
                             laplacian_smooth=0.1,
                             tolerance=0.1):
    '''
    Quadratic optimization of coil stream function according to minimal field energy,
    while keeping specified target field at target points within bounds.

    Optional Laplacian smoothing of inductance matrix.

    Parameters:
        meshobj: MeshWrapper object
        bfield_specification: list in which element is a dictionary containing a field specification
            each dict contains:
                C: Coupling matrix n_verts x n_verts
                target_field: n_r x 3
                abs_error: float
                rel_error: float
        objective: string or dict
            if string, either 'minimum_inductive_energy' or 'minimum_resistive_energy'
            if tuple, should contain: (a, b), where a and b are floats 0-1 describing the inductive and resitive weighting factors
        laplacian_smooth: float
        tolerance: float
    '''

    if objective == 'minimum_inductive_energy':
        objective = (1, 0)
    elif objective == 'minimum_resistive_energy':
        objective = (0, 1)

    #Initialize inequality constraint matrix and product
    constraint_matrix = np.zeros((0, len(meshobj.inner_verts)))
    constraint_product = np.zeros((0, ))

    #Populate inequality constraints with bfield specifications
    for spec in bfield_specification:

        #Limit C matrix to inner vertices (boundaries are kept at zero)
        inner_C = spec['C'][:, meshobj.inner_verts]

        #Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
        inner_C = inner_C.transpose((1, 0, 2))
        inner_C = inner_C.reshape((inner_C.shape[0], -1)).T

        #Apply relative error to bounds
        if spec['rel_error'] != 0:
            upper_bound = spec['target_field'] * (1 + np.sign(spec['target_field']) * spec['rel_error'])
            lower_bound = spec['target_field'] * (1 - np.sign(spec['target_field']) * spec['rel_error'])

        #Apply absolute error to bounds
        if spec['abs_error'] != 0:
            upper_bound = spec['target_field'] + spec['abs_error']
            lower_bound = spec['target_field'] - spec['abs_error']


        #Flatten to match C matrix
        upper_bound = upper_bound.flatten()
        lower_bound = lower_bound.flatten()


        #Stack upper and lower bounds into a single constraint
        stacked_bounds = np.hstack((lower_bound, -upper_bound))
        stacked_inner_C = np.vstack((-inner_C, inner_C))

        # Append specification to constraint matrix and product
        constraint_matrix = np.append(constraint_matrix, stacked_inner_C, axis=0)
        constraint_product = np.append(constraint_product, stacked_bounds, axis=0)


    #Limit L matrix to inner vertices
    inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

    #Limit R matrix to inner vertices
    inner_R = meshobj.resistance[meshobj.inner_verts][:, meshobj.inner_verts]


    #Linear part of QP problem not used, set to zero
    linear_part = np.zeros((len(meshobj.inner_verts), ))

    print('Scaling matrices before optimization. This requires singular value computation, hold on.')

    if laplacian_smooth != 0:

        #Limit Laplacian matrix to inner vertices (if used)
        inner_lapl = meshobj.laplacian.todense()[meshobj.inner_verts][:, meshobj.inner_verts]

        #Scale Laplacian matrix to same magnitude as L matrix and R matrix
        lapl_eigs = np.linalg.eigvalsh(-inner_lapl)

        L_eigs = np.linalg.eigvalsh(inner_L)

        R_eigs = np.linalg.eigvalsh(inner_R)

        scaled_lapl = np.max(np.abs(L_eigs)) / np.max(np.abs(lapl_eigs)) * -inner_lapl

        quadratic_term = (inner_L + laplacian_smooth * scaled_lapl)

    else:
        quadratic_term = inner_L

    #Scale whole quadratic term according to largest eigenvalue
    quad_eigs = np.linalg.eigvalsh(quadratic_term)
    quadratic_term /= np.max(np.abs(quad_eigs))


    #Compute, scale C matrix according to largest singular value
    u, s, vt = svds(constraint_matrix, k=1)

    #Also, scale constraints so max value is 1

    print('Solving quadratic programming problem using cvxopt...')
    I_inner, sol = cvxopt_solve_qp(P=quadratic_term,
                                   q=linear_part,
                                   G=constraint_matrix / s[0],
                                   h=constraint_product / np.max(constraint_product),
                                   tolerance=tolerance)

    #Build final I vector with zeros on boundary elements, scale by same singular value
    I = np.zeros((meshobj.inductance.shape[0], ))
    I[meshobj.inner_verts] = I_inner / s[0]

    return I, sol
