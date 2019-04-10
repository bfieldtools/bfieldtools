#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 17:08:27 2019

@author: Rasmus Zetter
"""
from scipy.optimize import minimize, LinearConstraint
import numpy as np

def optimize_streamfunctions(meshobj, target_field, target_axis, objective='energy', target_error=0.05):
    '''
    Optimize surface stream functions according to objective function.

    Example:
        Energy objective function minimizes field energy (and coil inductance)
        Linear constraint keeps B-field at target points within tolerated margin of error (set by target_error)
        NB! Constraint only applies to one B-field axis, set by target_axis.
        Other axes are freely minimized by objective function.

    '''
    if objective is 'energy':
        objective = lambda I, L: np.dot(np.dot(I, L), I)


    init_guess = np.random.random(size=(len(meshobj.inner_verts), ))

    lb = target_field * (1 - np.sign(target_field)* target_error)
    ub = target_field * (1 + np.sign(target_field)* target_error)

#    error_constraint = LinearConstraint(meshobj.C[:, meshobj.inner_verts, target_axis], lb, ub)

    inner_C = meshobj.C[:, meshobj.inner_verts, target_axis]


    lb_con = {"type":"ineq", "fun":lambda x: inner_C.dot(x) - lb}
    ub_con = {"type":"ineq", "fun":lambda x: -inner_C.dot(x) + ub}

    cons = [lb_con, ub_con]

    options = dict(disp=True, maxiter=2000)

    res = minimize(objective, init_guess, args=(meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]),
                                              constraints=cons, options=options)

    return res


import quadprog

def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)   # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]

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

def quadprog_streamfunctions(meshobj, target_field, target_axis, target_error=0.01):

    lb = target_field * (1 - np.sign(target_field)* target_error)
    ub = target_field * (1 + np.sign(target_field)* target_error)


    inner_C = meshobj.C[:, meshobj.inner_verts, target_axis]

    # Stack to turn upper and lower bound constraints into single constraint
    stacked_bounds = np.hstack((lb, -ub))
    stacked_inner_C = np.vstack((-inner_C, inner_C))

    inner_L = meshobj.inductance[meshobj.inner_verts][:, meshobj.inner_verts]
    linear_part = np.zeros((len(meshobj.inner_verts), ))

#    sol = quadprog_solve_qp(P=inner_L, q=linear_part, G=stacked_inner_C, h=stacked_bounds)
    sol = cvxopt_solve_qp(P=inner_L, q=linear_part, G=stacked_inner_C, h=stacked_bounds)

    I = np.zeros((meshobj.inductance.shape[0], ))
    I[obj.inner_verts] = sol

    return I
