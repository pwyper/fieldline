#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:06:50 2022

@author: nkqr53
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile

#set the tracing resolution (1 more than in fieldline.cnt)
xres = 401
yres = 401

ff = FortranFile('./PFLS_vert/qmap.0210003','r')
data = ff.read_reals(dtype=float)
qmap = data.reshape(xres,yres)
data = 0.0


#set any values less than 2 to 2 (min value Q can take)
qmap1 = qmap.copy()
lev = 2
qmap1[qmap1 < lev] = lev
qmap = qmap1

ff = FortranFile('./PFLS_vert/bbends.0210003','r')
data1 = ff.read_reals(dtype=float)
bends = data1.reshape(3,2,xres,yres)
data1 = 0.0


ff = FortranFile('./PFLS_vert/flen.0210003','r')
data = ff.read_reals(dtype=float)
flen = data.reshape(xres,yres)
data = 0.0

ff = FortranFile('./PFLS_vert/twist.0210003','r')
data = ff.read_reals(dtype=float)
twist = data.reshape(xres,yres)
data = 0.0

ff = FortranFile('./PFLS_vert/bmapping.0210003','r')
data = ff.read_reals(dtype=float)
bmapping = data.reshape(3,2,xres,yres)
data = 0.0

#adjust quantities for plotting

#set maximum value for flen
flen1 = flen.copy()
lev = 3
flen1[flen1 > lev] = lev
flen = flen1

#set maximum value for twist
twist1 = twist.copy()
lev = 3
twist1[twist1 > lev] = lev
twist = twist1

#plot values at the field line starting positions
#note: theta is in degrees lattitude
y = bmapping[0,0,0,:]
x = bmapping[1,0,:,0]
xx,yy = np.meshgrid(x,y)


#plotting
fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('binary')
mycmap2 = plt.get_cmap('bwr')

levels = np.linspace(np.log10(2),9,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(np.log10(qmap)),levels,cmap=mycmap)

#ax.set_aspect(1)
ax.set_ylabel(r'$R$',fontsize=18)
ax.set_xlabel(r'$\theta$',fontsize=18)
plt.savefig('./vert_Q.png')
#plt.show()


fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('binary')
mycmap2 = plt.get_cmap('bwr')

levels = np.linspace(0,3,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(flen),levels,cmap=mycmap)

#ax.set_aspect(1)
ax.set_ylabel(r'$R$',fontsize=18)
ax.set_xlabel(r'$\theta$',fontsize=18)
plt.savefig('./vert_fl_length.png')
#plt.show()


fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('bwr')

levels = np.linspace(-3,3,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(twist),levels,cmap=mycmap)     

#ax.set_aspect(1)
ax.set_ylabel(r'$R$',fontsize=18)
ax.set_xlabel(r'$\theta$',fontsize=18)
plt.savefig('./vert_Tw.png')
#plt.show()