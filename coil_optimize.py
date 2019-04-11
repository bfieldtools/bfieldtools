#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 17:08:27 2019

@author: Rasmus Zetter
"""

import numpy as np
import cvxopt
from cvxopt import matrix


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])
    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],))

def optimize_streamfunctions(meshobj, target_field, target_axis, target_error=0.1, laplacian_smooth=0.1):

    lb = target_field * (1 - np.sign(target_field)* target_error)
    ub = target_field * (1 + np.sign(target_field)* target_error)


    inner_C = meshobj.C[:, meshobj.inner_verts, target_axis]

    # Stack to turn upper and lower bound constraints into single constraint
    stacked_bounds = np.hstack((lb, -ub))
    stacked_inner_C = np.vstack((-inner_C, inner_C))


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
                   h=stacked_bounds)

    I = np.zeros((meshobj.inductance.shape[0], ))
    I[meshobj.inner_verts] = sol

    return I
