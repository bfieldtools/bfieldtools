import numpy as np
import cvxopt
from cvxopt import matrix
from scipy.sparse.linalg import svds
from scipy.linalg import eigh as largest_eigh

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, tolerance=1e-7):
    '''
    Use cvxopt to minimize
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
        return None, sol
    return np.array(sol['x']).reshape((P.shape[1],)), sol


def optimize_streamfunctions(meshobj, bfield_specification,
                             objective='minimum_inductive_energy',
                             laplacian_smooth=0.1,
                             tolerance=0.1):
    '''
    Quadratic optimization of coil stream function according to minimal field energy,
    while keeping specified target field at target points within bounds.

    Optional Laplacian smoothing of inductance matrix.

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


    #Linear part of QP problem not used, set to zero
    linear_part = np.zeros((len(meshobj.inner_verts), ))


    if objective == (1, 0):
        #Limit L matrix to inner vertices
        inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

        quadratic_term = inner_L

    elif objective == (0, 1):
        #Limit R matrix to inner vertices
        inner_R = meshobj.resistance[meshobj.inner_verts][:, meshobj.inner_verts]

        quadratic_term = inner_R

    else:
        #Limit L matrix to inner vertices
        inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

        #Limit R matrix to inner vertices
        inner_R = meshobj.resistance[meshobj.inner_verts][:, meshobj.inner_verts]

        print('Scaling inductance and resistance matrices before optimization. This requires eigenvalue computation, hold on.')

        max_eval_L = largest_eigh(inner_L, eigvals=(inner_L.shape[0]-1, inner_L.shape[0]-1))[0][0]
        max_eval_R = largest_eigh(inner_R, eigvals=(inner_L.shape[0]-1, inner_L.shape[0]-1))[0][0]

        scaled_R = max_eval_L / max_eval_R * inner_R

        quadratic_term = (objective[0] * inner_L  + objective[1] * scaled_R)


    #Scale whole quadratic term according to largest eigenvalue
    max_eval_quad = largest_eigh(quadratic_term, eigvals=(quadratic_term.shape[0]-1, quadratic_term.shape[0]-1))[0][0]

    quadratic_term /= max_eval_quad


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