#!/usr/bin/env python 
import os,sys  
import point_group as pg
from mysub import pmatrix,lstr
import numpy as np
import numpy.linalg as la
from periodica import structure
import fort_kpt_tools as fkpt
from scipy import *

def tetra_ibz(nk,ptg,struct):
  if type(nk)!=type(()):
    if [s for s in nk if type(s)!=int]:print "nk must be a tuple of integers";sys.exit()
  if not pg.grp.keys().count(ptg):print "point group %s not recognized"%ptg
  # get the structure...
  if type(struct)==type(''): pos=structure(struct)
  else: pos=struct

  def N(i,j,k):return i+(nk[0]+1)*(j+(nk[1]+1)*k)

  vec=np.transpose(pos.rvec.copy())
  print vec

  # form the direct coordinate symmetry elements...
  # NOTE: you must invoke round method or you will get 0 when converting to an integer array...
  G={}
  for g in pg.grp[ptg].G.keys():G[g]=np.dot(la.inv(vec),np.dot(pg.grp[ptg].G[g],vec)).round(5).astype(int)



  # define a mesh...
  ki=[(i,j,k) for k in range(nk[2]+1) for j in range(nk[1]+1) for i in range(nk[0]+1)]

#  # define the irreducible kpoints...
#  kptr=[]
#  for k in ki:
#    tt=N(*k);mo=k
#    for g in sorted(G.keys()):
#      temp=tuple([int(ss) if ss>=0 else int(ss+nk[s]) for s,ss in enumerate(np.dot(G[g],k))])
#      if N(*temp)<tt:mo=temp;tt=N(*temp)
#    kptr.append(mo)
#  kibz=sorted(set(kptr))

  # obtain the irreducible kpoints...
  nn=sorted(G.keys()) 
  if nn[0]!='E':print "Identity element MUST be first";sys.exit()
  print nk,[G[s] for s in nn]
  kptr,gptr=fkpt.get_kptr(nk,[G[s] for s in nn],ki)
  kptr=[tuple(s) for s in kptr]
  kibz=sorted(set(kptr))
  print kibz


  
  # now construct tetrahedra...
  # first we must determine which diagonal of the submesh cell to use...
  verts=((0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1))
  diags=((0,7),(1,6),(4,3))
  d=sum(np.dot(vec,diags[0][0]-diags[0][1])**2)
  # find the smallest diagonal...
  mdiag=sorted([(sum(np.dot(vec,np.array(verts[ss[0]])-np.array(verts[ss[1]]))**2),s) for s,ss in enumerate(diags)])[0][-1]
  mdiag=(verts[diags[mdiag][0]],verts[diags[mdiag][1]])
  # get the primitive tetrahedra...
  temp=[s for s in verts if s!=mdiag[0] and s!=mdiag[1]]
  ptet=[]
  for i,ii in enumerate(temp):
    for j,jj in enumerate(temp[i:]):
      if sum((np.array(ii)-np.array(jj))**2)==1:ptet.append(sorted([mdiag[0],mdiag[1],ii,jj])) 

  # now find all the irreducuble tetrahedra...
  ## this is what i originally did in python...
  #tet=[]
  #for i in range(nk[0]):
  #  for j in range(nk[1]):
  #    for k in range(nk[2]):
  #      for t in ptet: tet.append(tuple(sorted([kibz.index(kptr[N(*(np.array(s)+np.array([i,j,k])))]) for s in t])))
  ntet=nk[0]*nk[1]*nk[2]*6
  tetkptr,tet=fkpt.get_itet(nk,ptet,kibz,kptr,ntet)
  print tetkptr
  tetkptr=[tuple(s) for s in tetkptr]
  print tetkptr
  itet=sorted(set(tetkptr))
  #tetmult=[tet.count(s) for s in itet]
  tetmult=fkpt.get_tetmult(tetkptr,itet)


  #print len(ki),len(kibz),len(ki)/float(len(kibz))
  


  return kibz,itet,tetmult,tetkptr,tet,gptr


if __name__=="__main__":

  pos='POSCAR'
  nk=(6,6,3) 
  kibz,itet,tetmult,tetkptr,tet,gptr=tetra_ibz(nk,'c1',pos)
  norb=10;nki=100;ntet=6*6*3*6;
  #tetkptr=zeros((ntet,4),dtype=int)
  tet_idx=zeros(nki,dtype=int);eigval=zeros((nki,norb))
  cpdos=zeros((nki,norb))
  cpdos=fkpt.tetra_pnf(norb,nki,ntet,tet_idx,tetkptr)
  print kibz[tetkptr[0][0]-1],kibz[tetkptr[0][1]-1],kibz[tetkptr[0][2]-1],kibz[tetkptr[0][3]-1]
  print kibz[tetkptr[1][0]-1],kibz[tetkptr[1][1]-1],kibz[tetkptr[1][2]-1],kibz[tetkptr[1][3]-1]
  #import tb_tools as tb
  #pos='hi\n1\n1 0 0\n0 1 0\n0 0 1\n2\nd\n0 0 0 s\n0 0 0 s2'
  #pos='poscar'


  #kibz,itet,tetmult=tetra_ibz("10 10 10",'oh',pos)
  ##ham=tb.hamiltonian(pos,"%s/hartree_fock/test_ksum/rham_5band.py"%os.environ["HOME"],pg='c1')
  #ham=tb.hamiltonian(pos,"%s/hartree_fock/test_ksum/rham_allstates.py"%os.environ["HOME"],pg='oh')
  ##ham.dos_fort(nk=90) 


  
  def dos_fort(nbin=400,out='dos_tetra.out'):
    import fort_tb_tools as ftb
    import tempfort as tfort
    kmesh=(10,10,10)
    kibz,itet,tetmult,tetkptr,tet,gptr=tetra_ibz("%s %s %s"%kmesh,'oh',pos)
    print "kpoint info done" #,len(itet),len(tet)
    cdos,pdos,eps,deps={},{},{},{} 
    tt=sorted(ham.hopping.keys())
    hh=[ham.hopping[s] for s in tt]
    ss=[ham.sym.rrep[i] for i in sorted(ham.sym.rrep.keys())]
    #for i in ham.sym.rrep.keys():print i; pmatrix(ham.sym.rrep[i])
    #sys.exit()
    for sp in ham.spin:
      if 0: emin,emax,temp1,temp2=ftb.tetra_totdos(kmesh,nbin,tt,hh,kibz,itet,tetmult,12.5,13.0)
      else: emin,emax,temp1,temp2=tfort.tetra_pdos(kmesh,nbin,tt,hh,kibz,tetkptr,tet,gptr,ss)
      pdos[sp]=temp1;del temp1
      cdos[sp]=temp2;del temp2
      deps[sp]=(emax-emin)/nbin
      eps[sp]=[emin+deps[sp]*i for i in range(nbin+1)]
    del tt,hh

    # print out the dos...
    for sp in ham.spin:
      OUT=open(("_%s."%sp).join(out.split(".") if out.count(".") else [out,''] ),'w')
      #for i,ii in enumerate(pdos[sp]):OUT.write("%.6f "%(eps[sp][i]-(ham.fermi if ham.fermi else 0))+" %.6f "%ii+'\n')
      for i,ii in enumerate(pdos[sp]):OUT.write("%.6f "%(eps[sp][i]-(ham.fermi if ham.fermi else 0))+lstr(ii.tolist())+'\n')
      OUT.close()
      OUT=open('c'+("_%s."%sp).join(out.split(".") if out.count(".") else [out,''] ),'w')
      #for i,ii in enumerate(pdos[sp]):OUT.write("%.6f "%(eps[sp][i]-(ham.fermi if ham.fermi else 0))+" %.6f "%ii+'\n')
      for i,ii in enumerate(cdos[sp]):OUT.write("%.6f "%(eps[sp][i]-(ham.fermi if ham.fermi else 0))+lstr(ii.tolist())+'\n')
      OUT.close()


  #dos_fort()
