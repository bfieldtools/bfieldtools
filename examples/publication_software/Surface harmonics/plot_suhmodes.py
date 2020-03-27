#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:37:17 2020

@author: Rasmus Zetter
"""
from bfieldtools.conductor import Conductor
#from bfieldtools.suhtools import SuhBasis
import numpy as np

c = Conductor(mesh_file='/l/bfieldtools/examples/publication_software/Surface harmonics/arbitrary_surface2.stl', process=True, basis_name='suh', N_suh=15, fix_normals=True)


T_x = 1.5*np.pi/2
T_z = -1.02*np.pi
rotmat = np.array([[np.cos(T_z), -np.sin(T_z), 0, 0],
                   [np.sin(T_z), np.cos(T_z), 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]]) @ np.array([[1, 0, 0, 0],
                   [0, np.cos(T_x), -np.sin(T_x), 0],
                   [0, np.sin(T_x), np.cos(T_x), 0],
                   [0, 0, 0, 1]])



c.mesh.apply_transform(rotmat)


from mayavi import mlab
from mayavi.api import Engine
e = Engine()
e.start()



f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(750, 600))

c.suh_basis.plot(Nfuncs=c.basis.shape[1], dist=0.1, colormap='RdBu', figure=f, ncolors=256)

f.scene.z_plus_view()

#
for n in range(c.basis.shape[1]):
    surface = e.scenes[0].children[n].children[0].children[0].children[0]
#    surface.enable_contours = True
#    surface.contour.filled_contours = True
#    surface.contour.number_of_contours = 17
    surface.actor.mapper.interpolate_scalars_before_mapping = True


#f.scene.camera.parallel_projection=1
#mlab.view(0,160)
f.scene.camera.zoom(1.5)
#f.scene.camera.roll(270)

mlab.savefig('/l/bfieldtools/examples/publication_software/Surface harmonics/surface_harmonics.png', figure=f, magnification=4)

mlab.close()

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(750, 600))

c.plot_mesh(representation='wireframe', figure=f)
c.plot_mesh(opacity=0.2, figure=f)
#
#f.scene.camera.parallel_projection=1
#mlab.view(0,160)
f.scene.camera.zoom(1.5)
f.scene.z_plus_view()
#f.scene.camera.roll(270)

mlab.savefig('/l/bfieldtools/examples/publication_software/Surface harmonics/suhmesh.png', figure=f, magnification=4)
