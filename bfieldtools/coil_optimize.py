import numpy as np
import cvxopt
from cvxopt import matrix
from scipy.sparse.linalg import svds


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, tolerance=1e-7):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])

    cvxopt.solvers.options['abstol']=tolerance
    cvxopt.solvers.options['feastol']=tolerance
    cvxopt.solvers.options['reltol']=tolerance
#    cvxopt.solvers.options['refinement']=7

    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],)), sol

def optimize_streamfunctions(meshobj, target_field, target_axis,
                             target_error={'on_axis':0.05, 'off_axis':0.05, 'stray':0.05},
                             laplacian_smooth=0.1,
                             tolerance=0.1):
    '''
    Quadratic optimization of coil stream function according to minimal field energy,
    while keeping specified target field at target points within bounds.

    Optional Laplacian smoothing of inductance matrix.

    '''

#
    # Set lower and upper bound for stray field, all three axes
    lb_stray = np.repeat((-np.abs(np.repeat(target_field[0],len(meshobj.strayC), axis=0)) * target_error['stray'])[:, None], 3, axis=1).flatten()
    ub_stray = np.repeat((np.abs(np.repeat(target_field[0],len(meshobj.strayC), axis=0)) * target_error['stray'])[:, None], 3, axis=1).flatten()
#
#    #Limit stray field C matrix to inner vertices
    inner_strayC = meshobj.strayC[:, meshobj.inner_verts]
#
#    #Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
    inner_strayC = inner_strayC.transpose((1, 0, 2))
    inner_strayC = inner_strayC.reshape((inner_strayC.shape[0], -1)).T

    # Set lower and upper bounds on target axis
    lb_on_axis = target_field * (1 - np.sign(target_field) * target_error['on_axis'])
    ub_on_axis = target_field * (1 + np.sign(target_field) * target_error['on_axis'])

    # Set bounds on non-target axes
    lb_off_axis = -np.abs(target_field) * target_error['off_axis']
    ub_off_axis = np.abs(target_field) * target_error['off_axis']


    #Collect on-axis and off-axis bounds
    lb = np.full((lb_on_axis.shape[0], 3), fill_value=np.nan)
    ub = np.full((lb_on_axis.shape[0], 3), fill_value=np.nan)

    for ax in range(3):
        if ax is target_axis:
            lb[:, ax] = lb_on_axis
            ub[:, ax] = ub_on_axis
        else:
            lb[:, ax] = lb_off_axis
            ub[:, ax] = ub_off_axis

    #Reshape so that values are x1, y1, z1, x2, y2, z2, etc.
    lb = lb.flatten()
    ub = ub.flatten()

    #Limit C matrix to inner vertices
    inner_C = meshobj.C[:, meshobj.inner_verts]

    #Reshape so that values on axis 1 are x1, y1, z1, x2, y2, z2, etc.
    inner_C = inner_C.transpose((1, 0, 2))
    inner_C = inner_C.reshape((inner_C.shape[0], -1)).T

    # Stack to turn upper and lower bound constraints into single constraint
    stacked_bounds = np.hstack((lb, -ub, lb_stray, -ub_stray))
    stacked_inner_C = np.vstack((-inner_C, inner_C, -inner_strayC, inner_strayC))

    #Limit L matrix to inner vertices
    inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]

    #Linear part of QP problem not used, set to zero
    linear_part = np.zeros((len(meshobj.inner_verts), ))



    print('Scaling matrices before optimization. This requires singular value computation, hold on.')

    if laplacian_smooth != 0:

        #Limit Laplacian matrix to inner vertices (if used)
        inner_lapl = meshobj.laplacian.todense()[meshobj.inner_verts][:, meshobj.inner_verts]

        #Scale Laplacian matrix to same magnitude as inductance i.e. L
        lapl_eigs = np.linalg.eigvalsh(-inner_lapl)
        L_eigs = np.linalg.eigvalsh(inner_L)

        scaled_lapl = np.max(np.abs(L_eigs))/np.max(np.abs(lapl_eigs))*-inner_lapl

        quadratic_term = (inner_L + laplacian_smooth * scaled_lapl)

    else:
        quadratic_term = inner_L

    #Scale whole quadratic term according to largest eigenvalue
    quad_eigs = np.linalg.eigvalsh(quadratic_term)
    quadratic_term /= np.max(np.abs(quad_eigs))


    #Compute, scale C matrix according to largest singular value
    u, s, vt = svds(stacked_inner_C, k=1)

    #Also, scale constraints so max value is 1

    print('Solving quadratic programming problem using cvxopt...')
    I_inner, sol = cvxopt_solve_qp(P=quadratic_term,
                   q=linear_part,
                   G=stacked_inner_C/s[0],
                   h=stacked_bounds/np.max(target_field),
                   tolerance=tolerance)

    #Build final I vector with zeros on boundary elements, scale by same singular value
    I = np.zeros((meshobj.inductance.shape[0], ))
    I[meshobj.inner_verts] = I_inner / s[0]

    return I, sol
