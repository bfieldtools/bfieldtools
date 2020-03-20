

from bfieldtools import mesh_class
from bfieldtools.utils import load_example_mesh
import pkg_resources

import pytest

import numpy as np
from numpy.testing import (assert_array_almost_equal, assert_array_equal,
                           assert_allclose, assert_equal)



def _fake_conductor(mesh_name='10x10_plane', **kwargs):
    """
    Creates an example 'fake' Conductor object
    """
    return mesh_class.Conductor(mesh_obj=load_example_mesh(mesh_name), **kwargs)


def _fake_streamfunction(vals=None, **kwargs):
    """
    Creates an example 'fake' StreamFunction object
    """
    if vals:
        return mesh_class.StreamFunction(vals, _fake_conductor())
    else:
        
        return mesh_class.StreamFunction()


def test_conductor_creation():
    """
    Tests different ways to create Conductor objects
    """
    for test_mesh in ['unit_sphere', 'unit_plane', 'plane_with_holes']:
        for basis_name in ['suh', 'inner', 'vertex']:
            c = _fake_conductor(mesh_name=test_mesh, basis_name=basis_name)
            

def test_conductor_basis_change():
    
    c = _fake_conductor(basis_name='vertex')
    
    assert c
    
    c.set_basis('inner')
    
    assert c
    
    c.set_basis('suh')
    
    assert c
    
    
def test_streamfunction_creation():
    """
    Test creating StreamFunctions with different bases
    """
    
    mesh_class.StreamFunction(np.zeros((10,)), _fake_conductor(basis_name='suh', N_suh=10))
    mesh_class.StreamFunction(np.zeros((584,)), _fake_conductor(basis_name='inner'))
    mesh_class.StreamFunction(np.zeros((676,)), _fake_conductor(basis_name='vertex'))
    
    
    
def test_conductor_attributes():
    """
    tests Conductor attributes
    """
    pass
    

def test_streamfunction_attributes():
    """
    """
    pass


    