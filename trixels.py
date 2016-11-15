#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Python example of how to use trixels.x

This small python scrip makes a comparison between FFT and the trixel method
of Fourier transformation presented by Brinch & Dullemond, 2014, mnras, 440, 3285

This model used in this example is a thin ring.

"""

__author__     = "Christian Brinch"
__copyright__  = "Copyright 2014-2016"
__credits__    = ["Christian Brinch"]
__license__    = "AFL 3.0"
__version__    = "1.0"
__maintainer__ = "Christian Brinch"
__email__      = "brinch@nbi.ku.dk"

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.delaunay as triang
from matplotlib.backends.backend_pdf import PdfPages
from subprocess import Popen, PIPE

def radius(i,j):
    return (float(i)-npix/2.)**2 + (float(j)-npix/2.)**2

def lloyd(xn,yn):
    f = open('points.in','w')
    f.write("2 rbox %d D2\n" % len(xn))
    f.write("%d\n" % len(xn))
    for p in range(len(xn)):
        f.write("%e " % xn[p])
        f.write("%e\n" % yn[p])
    f.close()

    for count in range(20):
        cmd = 'cat points.in | qvoronoi s o QJ'
        sub_process = Popen(cmd, shell=True,stdout=PIPE,stderr=PIPE)
        output = sub_process.communicate()
        lines = [line.split() for line in output[0].split('\n') if line]

        setsize = len(lines)
        nv = int(lines[1][0])
        coord = lines[2:2+nv]
        ref = lines[2+nv:]

        for P in range(len(xn)):
            cx = 0.
            cy = 0.

        for V in range(int(ref[P][0])):
            cx += float(coord[int(ref[P][V+1])][0])
            cy += float(coord[int(ref[P][V+1])][1])

        n = int(ref[P][0])
        xn[P] = xn[P] - (xn[P]-cx/n)/15.
        yn[P] = yn[P] - (yn[P]-cy/n)/15.

        f = open('points.in','w')
        f.write("2 rbox %d D2\n" % len(xn))
        f.write("%d\n" % len(xn))
        for p in range(len(xn)):
            f.write("%e " % xn[p])
            f.write("%e\n" % yn[p])
        f.close()

    return triang.delaunay(xn,yn)




pp = PdfPages('transforms.pdf')
fig = plt.figure(1)

ss = 8
immin = -1
immax = 3.5
npix = 51
rout = 72
rin = 36

image = np.array([ [10. if radius(i,j) < rout and radius(i,j) > rin else 0.0001 \
         for i in range(npix) ] for j in range(npix) ])


# Contour pixel image
x = np.linspace(0, npix, npix*10)
X,Y = np.meshgrid(x,x)
ax=plt.subplot(221)
ax.set_xlabel('x')
ax.set_ylabel('y')
image_big = np.array([[image[i/10,j/10] for i in range(10*npix)] \
                     for j in range(10*npix)])
ax = plt.contourf(X, Y, image_big, levels=[0,1,2,3,4,5,6,7,8,9,10], \
                  cmap=plt.cm.Reds)
ax = plt.colorbar()
for i in range(npix):
  ax = plt.plot([i,i],[0,npix],color='grey',lw=0.1)
for i in range(npix):
  ax = plt.plot([0,npix],[i,i],color='grey',lw=0.1)



# Fast Fourier transform pixel image and contour plots
FT = np.fft.fft2(image)
n = FT.shape[0]
freq = np.fft.fftshift(np.fft.fftfreq(n,1))
real = np.max(FT.real)
imag = np.max(FT.imag)


ax=fig.add_subplot(222)
ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_xlim(-0.5,0.5)
ax.set_ylim(-0.5,0.5)
FTshift = np.abs(np.fft.fftshift(np.sqrt(FT**2)))
im = ax.imshow(np.log10(FTshift), interpolation='nearest', origin='lower',\
               extent=[np.min(freq),np.max(freq),np.min(freq),np.max(freq)],\
               vmin=immin, vmax=immax)
ax = plt.colorbar(im)





# Make a weighted random trixel grid
x = [[],[],[]]
for i in range(npix**2):
  flag=True
  while(flag):
    tx = np.random.uniform(0,npix)
    ty = np.random.uniform(0,npix)

    if radius(ty,tx) < rout and radius(ty,tx) > rin:
        tz = 10.
    else:
        tz = 0.

    if radius(ty,tx) < (rout+0.2*rout) and radius(ty,tx) > (rin-0.2*rin):
        flag = False
    else:
        if(np.random.uniform(0,1,1) < 1.e-3):
            flag = False

  x[0].append(tx)
  x[1].append(ty)
  x[2].append(tz)


try:
    cens,edg,tri,neig = lloyd(x[0],x[1])
except:
    print "Qhull not found. Skipping lloyd smoothing."
    cens,edg,tri,neig = triang.delaunay(x[0],x[1])


ax = fig.add_subplot(223)
ax.set_xlim(0,npix)
ax.set_ylim(0,npix)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax = plt.tripcolor(x[0], x[1], tri, x[2], shading='gouraud', cmap=plt.cm.Reds)
ax = plt.colorbar()
ax = plt.tripcolor(x[0], x[1], tri, x[2], shading='flat',alpha=0.5, \
     edgecolors='black', lw=0.1, cmap=plt.cm.Reds)


FILE = open("temp_out","w")
FILE.write(str(len(x[0]))+"\n")
FILE.write(str(tri.size/3)+"\n")
FILE.write(str(npix*ss)+"\n")
for i in range(len(x[0])):
      FILE.write(str(x[0][i])+"\n")
for i in range(len(x[0])):
      FILE.write(str(x[1][i])+"\n")
for i in range(len(x[0])):
      FILE.write(str(x[2][i])+"\n")
for i in range(tri.size/3):
  for j in range(3):
      FILE.write(str(tri[i,j])+"\n")
FILE.close()

os.system('./trixels.x')
vis = np.loadtxt("temp_in")
coords = np.loadtxt("temp_co")
os.system('rm -rf temp_out')
os.system('rm -rf temp_in')
os.system('rm -rf temp_co')

ndim=npix*ss
vis.shape=(ndim,ndim,2)
vis[:,:,0]=np.transpose(vis[:,:,0])
vis[:,:,1]=np.transpose(vis[:,:,1])

ax=plt.subplot(224)
ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_xlim(-0.5,0.5)
ax.set_ylim(-0.5,0.5)
freq = np.linspace(-(npix/2.-1)/npix, (npix/2.-1)/npix, num=ndim)
FTshift = np.sqrt(vis[:,:,0]**2+vis[:,:,1]**2)
vis=vis*np.pi/2.
im=ax.imshow(np.log10(FTshift), interpolation='nearest', origin='lower', \
             extent=[np.min(freq),np.max(freq),np.min(freq),np.max(freq)],\
             vmin=immin, vmax=immax)
ax=plt.colorbar(im)


plt.savefig(pp, format='pdf')
pp.close()
