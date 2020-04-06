# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:52:32 2020

@author: Rasmus Zetter
"""
import numpy as np
from trimesh.path import Path3D
from trimesh.path.entities import Line

from .viz import plot_3d_current_loops
from . import line_magnetics



class LinePath(Path3D):
    '''
    Class that inherits Trimesh.path.Path3D for handling discretized current loops

    '''

    def __init__(self, loops):
        '''
        Init with inheritance 
        '''
        vertices = np.zeros((0, 3))
        entities = []
        
        for loop in loops:
            
            points = np.append(np.arange(0, len(loop)), 0) + len(vertices)
            entities.append(Line(points))
            
            vertices = np.append(vertices, loop, axis=0)      
        
        
        Path3D.__init__(self, entities, vertices)
        self.attr3 = 'three'
        
    def plot_loops(self, **kwargs):
        '''
        Plots loops in 3D using mayavi, see viz.plot_3d_current_loops for more details
        
        Parameters
        ----------
        colors: str
        
        
        '''    

        if 'tube_radius' not in kwargs:
            kwargs['tube_radius'] = 0.005*np.linalg.norm(self.bounds)

        figure = plot_3d_current_loops([e.discrete(self.vertices) for e in self.entities], 
                                       **kwargs)
        return figure
    
    
    def magnetic_field(self, points, separate_loops=False):
        '''
        Compute magnetic field in some point due to a unit current in the loops
        
        Parameters
        ----------
        points: array (N_p, 3)
        
        separate_loops: Boolean
            	If True, don't combine contributions of separate loops
        
        Returns
        -------
        
        Bfield: array (N_p, 3) or array (N_loops, N_p, 3)
        
        '''
        Bfield = np.zeros((len(self.entities), len(points), 3))
        for ii, loop in enumerate(self.entities):
            Bfield[ii] = line_magnetics.magnetic_field(self.vertices[loop.points], points)
            
        if not separate_loops:
            Bfield = np.sum(Bfield, axis=0)
            
        return Bfield
        
    
        
    def vector_potential(self, points, separate_loops=False, **kwargs):
        '''
        Compute magnetic vector potential in some point due to a unit current in the loops
        
        Parameters
        ----------
        points: array (N_p, 3)
        
        separate_loops: Boolean
            	If True, don't combine contributions of separate loops
        
        Returns
        -------
        
        Aield: array (N_p, 3) or array (N_loops, N_p, 3)
        '''
        
        Afield = np.zeros((len(self.entities), len(points), 3))
        for ii, loop in enumerate(self.entities):
            Afield[ii] = line_magnetics.vector_potential(self.vertices[loop.points[:-1]], points, **kwargs)[0, :, :]
            
        if not separate_loops:
            Afield = np.sum(Afield, axis=0)
            
        return Afield
    
    
    def scalar_potential(self, points, separate_loops=False, **kwargs):
        '''
        Compute magnetic scalar potential in some point due to a unit current in the loops
        
        Parameters
        ----------
        points: array (N_p, 3)
        
        separate_loops: Boolean
            	If True, don't combine contributions of separate loops
        
        Returns
        -------
        
        Ufield: array (N_p,) or array (N_loops, N_p)
        '''

        Ufield = np.zeros((len(self.entities), len(points)))
        for ii, loop in enumerate(self.entities):
            Ufield[ii] = line_magnetics.scalar_potential(self.vertices[loop.points], points, **kwargs)
            
        if not separate_loops:
            Ufield = np.sum(Ufield, axis=0)
            
        return Ufield      
        
        
        