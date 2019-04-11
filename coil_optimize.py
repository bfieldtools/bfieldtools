#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 17:08:27 2019

@author: Rasmus Zetter
"""

import numpy as np
import cvxopt
from cvxopt import matrix


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, tolerance=1e-7):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])

    cvxopt.solvers.options['abstol']=tolerance

    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],))

def optimize_streamfunctions(meshobj, target_field, target_axis,
                             target_error=0.1,
                             laplacian_smooth=0.1,
                             tolerance=1e-7):
    '''
    Quadratic optimization of coil stream function according to minimal field energy,
    while keeping specified target field at target points within bounds.

    Optional Laplacian smoothing of inductance matrix.

    '''

    # Set lower and upper bounds on target axis
    lb_on_axis = target_field * (1 - np.sign(target_field) * target_error['on_axis'])
    ub_on_axis = target_field * (1 + np.sign(target_field) * target_error['on_axis'])

    # Set bounds on non-target axes
    lb_off_axis = -np.abs(target_field) * target_error['off_axis']
    ub_off_axis = np.abs(target_field) * target_error['off_axis']

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
    stacked_bounds = np.hstack((lb, -ub))
    stacked_inner_C = np.vstack((-inner_C, inner_C))

    #Limit L matrix to inner vertices
    inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]


    linear_part = np.zeros((len(meshobj.inner_verts), ))


    inner_lapl = meshobj.laplacian.todense()[meshobj.inner_verts][:, meshobj.inner_verts]

    #Scale laplacian matrix to same magnitude as inductance
    lapl_eigs = np.linalg.eigvalsh(-inner_lapl)
    ind_eigs = np.linalg.eigvalsh(inner_L)

    scaled_lapl = np.max(np.abs(ind_eigs))/np.max(np.abs(lapl_eigs))*-inner_lapl



    sol = cvxopt_solve_qp(P=inner_L + laplacian_smooth * scaled_lapl,
                   q=linear_part,
                   G=stacked_inner_C,
                   h=stacked_bounds,
                   tolerance=tolerance)

    I = np.zeros((meshobj.inductance.shape[0], ))
    I[meshobj.inner_verts] = sol

    return I
