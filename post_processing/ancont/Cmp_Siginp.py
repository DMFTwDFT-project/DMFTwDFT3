#!/usr/bin/env python
from scipy import *

Sig=loadtxt('sig.inpx')
norb=(shape(Sig)[1]-1)/2
SigMdc=zeros(norb)
for i in range(norb): SigMdc[i]=Sig[-1][1+2*i]
savetxt('Sigoo.out',[SigMdc[:norb]])
for i in range(norb): Sig[:,1+2*i]-=SigMdc[i]
fline=open('sig.inpx','r').readline()[:-1]
print fline
savetxt('sig.inp',Sig[:,:1+2*norb],header=fline)
