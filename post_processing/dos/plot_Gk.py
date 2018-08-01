#!/usr/bin/env python

from scipy import *
import matplotlib
matplotlib.use('ps')
from matplotlib.font_manager import fontManager, FontProperties
from pylab import *
import scipy.interpolate
import re

vmm = [0,5.0]
nk=0
SKP=[]
SKPoints=[]
distk=[]
kpts=[]
fi=open('klist.dat','r')
for line in fi.readlines():
   line=line.split()
   distk.append(float(line[0]))
   kpts.append([float(line[1]),float(line[2]),float(line[3])])
   if len(line)==5:
      SKP.append(float(line[0]))
      SKPoints.append(line[4])
#print distk
#print kpts
#print SKP
#print SKPoints
fi.close() 
fi=open('ksum.input','r')
numk=int(fi.readline())
nom=int(fi.readline())
fi.close()
A_k=[]
dist_k=[]
om=[]
kpts=[]
fi=open('Gk.out','r')
for i in range(numk):
   kpts.append(map(float,fi.readline().split()[1:]))
   A_k.append([]) 
   om.append([]) 
   for j in range(nom):
      line=map(float,fi.readline().split())
      A_k[i].append(-1*line[2]/3.14159265)
      om[i].append(line[0])
   A_k[i]=A_k[i][::-1]
fi.close()
A_k=array(matrix(array(A_k)[::-1]).T)
#print A_k
#exit()

#om=[]
#for line in file.readlines()[::-1]:
#   dat = map(float,line.split())
#   om.append(dat[0])
#   t_A_k.append(dat[1])
#A_k.append(t_A_k)
#
#for i in range(1,num_k+1):
#   file=open('G_k.out.%03d'%(i),'r')
#   t_A_k=[]
#   om=[]
#   for line in file.readlines()[::-1]:
#      dat = map(float,line.split())
#      om.append(dat[0])
#      t_A_k.append(dat[1])
#   A_k.append(t_A_k)
#A_k=array(matrix(array(A_k)).T)
# 
#k_eig=[]
#eigval0=[]
#file=open('eigval0.dat','r')
#for i,line in enumerate(file.readlines()):
#   dat = map(float,line.split())
#   eigval0.append([])
#   k_eig.append(dat[0])
#   for j in range(len(dat)-1):
#      eigval0[i].append(dat[j+1])
#eigval0=array(eigval0)
#k_eig=2*array(k_eig)/(float(k_eig[-1])+0.0)


(ymin,ymax) = (om[0][0],om[0][-1])
#print ymin,ymax
(xmin,xmax) = (distk[0],distk[-1])#(0, (shape(A_k)[1]-1)/(numk+0.0)*3.14159265)
#print shape(A_k),xmin,xmax
#vmm = [0,max(map(max,A_k))*itensity]

im=imshow(A_k,cmap=cm.hot, vmin=vmm[0], vmax=vmm[1], extent=[xmin,xmax,ymin,ymax],aspect='auto')
#im=imshow(A_k,cmap=cm.hot, interpolation='bilinear', vmin=vmm[0], vmax=vmm[1], extent=[0.0,10.0,-1.0,1.0])
colorbar(im,orientation='vertical',pad=0.05,shrink=1.0,ticks=arange(0,10.0,1.0))
#ax1=subplot(111)
#for i in range(18,shape(eigval0)[1]):
#   ax1.plot(k_eig,eigval0[:,i]-7.22220722842 ,color='green')
##ax1.plot([k_eig[319],k_eig[319]], [ymin,ymax], 'w-')
#xticks([0,1,2],['$\Gamma$','X','M'])
xticks(SKP,SKPoints)

axhline(y=0,color='black',ls='--')
xlabel('k-path',fontsize='xx-large')
ylabel('Energy',fontsize='xx-large')
#ax1.set_ylim(-1,0.5)


show()
savefig('A_k.eps')
exit()


