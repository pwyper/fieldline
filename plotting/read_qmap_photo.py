#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:06:50 2022

@author: nkqr53
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile


xres = 401
yres = 401

ff = FortranFile('./PFLS_photo/qmap.0210003','r')
data = ff.read_reals(dtype=float)
qmap = data.reshape(xres,yres)
data = 0.0

#set any values less than 2 to 2 (min value Q can take)
qmap1 = qmap.copy()
lev = 2
qmap1[qmap1 < lev] = lev
qmap = qmap1

ff = FortranFile('./PFLS_photo/bbends.0210003','r')
data1 = ff.read_reals(dtype=float)
bends = data1.reshape(2,3,xres,yres)
data1 = 0.0


ff = FortranFile('./PFLS_photo/flen.0210003','r')
data = ff.read_reals(dtype=float)
flen = data.reshape(xres,yres)
data = 0.0

ff = FortranFile('./PFLS_photo/twist.0210003','r')
data = ff.read_reals(dtype=float)
twist = data.reshape(xres,yres)
data = 0.0

ff = FortranFile('./PFLS_photo/bmapping.0210003','r')
data = ff.read_reals(dtype=float)
bmapping = data.reshape(3,2,xres,yres)
data = 0.0

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

x = bmapping[1,0,0,:]
y = bmapping[2,0,:,0]

#x = np.linspace(-30,30,xres)
#y = np.linspace(-30,30,yres)
xx,yy = np.meshgrid(x,y)




fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('binary')
mycmap2 = plt.get_cmap('bwr')

levels = np.linspace(np.log10(2),12,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(np.log(qmap)),levels,cmap=mycmap)
ax.contour(xx,yy,np.transpose(bends[0,0,:,:]),[0],linewidths=1,linestyles='dashed',colors='black')
ax.contourf(xx,yy,np.transpose(bends[0,0,:,:]),levels2,cmap=mycmap2,alpha=0.1)


ax.set_aspect(1)
ax.set_ylabel(r'$\theta$',fontsize=18)
ax.set_xlabel(r'$\phi$',fontsize=18)
plt.savefig('./photo_Q.png')
#plt.show()


fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('binary')
mycmap2 = plt.get_cmap('bwr')

levels = np.linspace(0,3,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(flen),levels,cmap=mycmap)
ax.contour(xx,yy,np.transpose(bends[0,0,:,:]),[0],linewidths=1,linestyles='dashed',colors='black')
#ax.contourf(xx,yy,np.transpose(bends[0,0,:,:]),levels2,cmap=mycmap2,alpha=0.1)


ax.set_aspect(1)
ax.set_ylabel(r'$\theta$',fontsize=18)
ax.set_xlabel(r'$\phi$',fontsize=18)

plt.savefig('./photo_fl_length.png')
#plt.show()


fig, ax = plt.subplots(constrained_layout=False,figsize=(6,6))

mycmap = plt.get_cmap('bwr')

levels = np.linspace(-3,3,50)
levels2 = np.linspace(-40,40,50)

ax.contourf(xx,yy,np.transpose(twist),levels,cmap=mycmap)
ax.contour(xx,yy,np.transpose(bends[0,0,:,:]),[0],linewidths=1,linestyles='dashed',colors='black')
#ax.contourf(xx,yy,np.transpose(bends[0,0,:,:]),levels2,cmap=mycmap2,alpha=0.1)


ax.set_aspect(1)
ax.set_ylabel(r'$\theta$',fontsize=18)
ax.set_xlabel(r'$\phi$',fontsize=18)

plt.savefig('./photo_Tw.png')
#plt.show()