#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 13:39:28 2019

@author: Rasmus Zetter
"""

import numpy as np

class CouplingMatrix:
    '''
    General-use class that contains a data array (a coupling matrix)
    and a bookkeeping list of computed points.

    When called, returns the coupling matrix for queried points.
    If some output has already been computed, use pre-computed values instead
    and only compute missing parts.

    TODO: What if parent is modified after being passed? How does this pass-by-reference act?

    '''
    def __init__(self, parent, function):

        #Bookkeeping array, which points are already computed
        #Indexed in same order as the matrix
        self.points = np.array([])

        self.matrix = np.array([])

        self.parent = parent
        self.function = function


    def __call__(self, points, *fun_args):
        '''
        Returns the output of self.function(self.parent, points).
        If some output has already been computed, use pre-computed values instead
        and only compute missing parts.

        Parameters
        ----------
            points: (N_points, ... ) numpy array
                Array containing query points
            *fun_args: additional arguments
                Optional, additional arguments that are passed to self.function

        Returns
        -------
            (N_points, ...) numpy array
        '''

        if len(self.points) == 0:
            self.matrix = self.function(self.parent.mesh, points, *fun_args)
            self.points = points

            return self.matrix

        else:
            #Check which points exist and which need to be computed
            m_existing_point_idx, p_existing_point_idx = np.where((self.points == points).all(axis=1))


            #Check if there are missing points
            missing_point_idx = np.arange(0, len(points))[~p_existing_point_idx]

            #If there are missing points, compute and add
            if missing_point_idx:
                missing_points = points[missing_point_idx]

                new_matrix_elems = self.function(self.parent.mesh, missing_points, *fun_args)


                #Append newly computed point to coupling matrix, update bookkeeping

                self.points.append(missing_points)
                self.matrix = np.vstack((self.matrix, new_matrix_elems))

                #Re-compute indices of queried points, now that all should exist
                m_existing_point_idx, p_existing_point_idx = np.where((self.points == points).all(axis=1))


            return self.matrix[m_existing_point_idx]


