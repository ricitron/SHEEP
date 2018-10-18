#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 09:43:10 2018

@author: max
"""

mesh_file = '/home/max/Downloads/mesh-00000.h5'
solution_file = '/home/max/Downloads/solution-00174.h5'

import h5py as h5
import numpy as np

## Load the mesh
mesh = h5.File(mesh_file,'r')
nodes = np.array(mesh['nodes'])

solution = h5.File(solution_file,'r')
basal_layer = np.array(solution['basal_layer'])
T = np.array(solution['T'])

## identify duplicate nodes - note that this is a bit sloppy because it presumes that the nodal values are also duplicated, which may not be true for DG fields
unique_nodes, unique_indices, unique_inverse, unique_counts = np.unique(nodes,axis=0,return_index=True,return_inverse=True,return_counts=True)

nx = np.unique(nodes[:,0]).size
ny = np.unique(nodes[:,1]).size

x_new = np.reshape(unique_nodes[:,0],[nx,ny])
y_new = np.reshape(unique_nodes[:,1],[nx,ny])

T_new = np.reshape(T[unique_indices],[nx,ny])
composition_new = np.reshape(basal_layer[unique_indices],[nx,ny])

import matplotlib.pyplot as plt
plt.figure()
plt.pcolor(x_new,y_new,T_new )
plt.colorbar()
plt.axis('equal')
plt.show()

plt.figure()
plt.pcolor(x_new,y_new,composition_new )
plt.colorbar()
plt.axis('equal')
plt.show()

## less sloppy approach
import scipy.spatial as sps
tree = sps.cKDTree(nodes)
# note - k=4 is appropriate for 2D only. Choice of 1 for distance upper bound is appropriate for lengths in m.
d,idx = tree.query(nodes,k=4,distance_upper_bound=1.0)
# pull out just the rows of idx corresponding to the unique indices
idx = idx[unique_indices,:]
d = d[unique_indices,:]

composition_new2 = np.zeros_like(unique_indices,dtype=float)
for i in range(idx.shape[0]):
    mask = d[i,:] <= 1.
    composition_new2[i] = np.mean(basal_layer[idx[i,mask]])
composition_new2 = np.reshape(composition_new2,[nx,ny])

plt.figure()
plt.pcolor(x_new,y_new,composition_new2 )
plt.colorbar()
plt.axis('equal')
plt.show()