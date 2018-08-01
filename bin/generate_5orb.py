import os,sys
from scipy import *
import scipy.linalg as lin
import copy

def Construct_basis(N_states):
   ####### Construct Fock basis ###########
   occ=[]
   for i in range(N_states):
      occ.append([0,1])
   basis_big=[]
   for i in range(N_states):
      basis_big.append([])
      for j in occ[i]:
         if i==0: basis_big[i].append([j])
         else:
            for k in basis_big[i-1]:
               basis_big[i].append(k+[j])

   basis=copy.deepcopy(basis_big[-1])
   return basis

def Construct_atomic_energy(J,e_d,norb,nspin,basis):
   E_atom=[];n_atom=[];Sz=[]
   for i,bas in enumerate(basis):
      E_atom.append([]); n_atom.append([]); Sz.append([])
      for k in range(len(bas)):
         if len(bas[k]) != norb*nspin or len(bas[k]) != len(e_d): print "Something is wrong!"; exit()
         N_up=0; N_dn=0; N_updn=0; E_orb=0.0
         bs=bas[k]
         for iorb in range(norb):
            N_up+=bs[2*iorb]; N_dn+=bs[2*iorb+1]; N_updn+=bs[2*iorb]*bs[2*iorb+1]
            E_orb+=(e_d[2*iorb])*bs[2*iorb]+(e_d[2*iorb+1])*bs[2*iorb+1]
         E_atom[i].append(E_orb-J*(N_up*N_dn+N_dn*N_up-2*N_updn)-1.5*J*(N_up*(N_up-1)+N_dn*(N_dn-1)))
         Sz[i].append(0.5*(N_up-N_dn))
         n_atom[i].append(N_up+N_dn)
   return (E_atom,n_atom,Sz)

def Construct_atomic_energy2(Ueff,e_d,norb,nspin,basis):
   E_atom=[];n_atom=[];Sz=[]
   for i,bas in enumerate(basis):
      E_atom.append([]); n_atom.append([]); Sz.append([])
      for k in range(len(bas)):
         if len(bas[k]) != norb*nspin or len(bas[k]) != len(e_d): print "Something is wrong!"; exit()
         N_up=0; N_dn=0; N_updn=0; E_orb=0.0
         bs=bas[k]
         for iorb in range(norb):
            N_up+=bs[2*iorb]; N_dn+=bs[2*iorb+1]; N_updn+=bs[2*iorb]*bs[2*iorb+1]
            E_orb+=(e_d[2*iorb])*bs[2*iorb]+(e_d[2*iorb+1])*bs[2*iorb+1]
            for jorb in range(norb):
               E_orb+=0.5*(Ueff[0][iorb,jorb]*bs[2*iorb]*bs[2*jorb]+Ueff[1][iorb,jorb]*bs[2*iorb]*bs[2*jorb+1])
               E_orb+=0.5*(Ueff[1][iorb,jorb]*bs[2*iorb+1]*bs[2*jorb]+Ueff[0][iorb,jorb]*bs[2*iorb+1]*bs[2*jorb+1])
         E_atom[i].append(E_orb)
         Sz[i].append(0.5*(N_up-N_dn))
         n_atom[i].append(N_up+N_dn)
   return (E_atom,n_atom,Sz)

def Off_J(basis,basis_pair,E_atom):
#   diag_E_atom=copy.deepcopy(E_atom)
   diag_E_atom=[]
   eigv=[]
   H_imp=zeros((2,2),dtype=float)
   for i,bs in enumerate(basis):
      if len(bs)==2:
         for bs2 in basis_pair:
            chk=0
            if bs[0]==bs2[0]: H_imp[0,0]=E_atom[i][0]; H_imp[1,1]=E_atom[i][1]; H_imp[1,0]=bs2[2]; H_imp[0,1]=bs2[2]; chk=1; break 
         if chk==0: print "Basis_pair is wrong!"; exit()
         (w,vr)=lin.eigh(H_imp)
         diag_E_atom.append([w[0],w[1]])
    #     if vr[0,0]<0: eigv.append([[-vr[0,0],-vr[1,0]],[vr[0,1],vr[1,1]]])
    #     elif vr[0,1]<0: eigv.append([[vr[0,0],vr[1,0]],[-vr[0,1],-vr[1,1]]])
         eigv.append([[vr[0,0],vr[1,0]],[vr[0,1],vr[1,1]]])
      else: 
         diag_E_atom.append(E_atom[i])
         eigv.append([[1.0]])
   return (diag_E_atom,eigv)

def Construct_Operator(basis,eigv,Op):
   Op_val=[]; 
   for i in range(len(Op)): Op_val.append([])
   for i in range(len(Op)):
      for j,bs1 in enumerate(basis):
         #print bs1
         if len(bs1)==1:
            Op_val[i].append([eigv[j][0][0]**2*sum(bs1[0][ii] for ii in Op[i])])
         elif len(bs1)==2:
            Op_val[i].append([sum(eigv[j][0][jj]*eigv[j][0][jj]*sum(bs1[jj][ii] for ii in Op[i]) for jj in range(2)),\
                              sum(eigv[j][0][jj]*eigv[j][1][jj]*sum(bs1[jj][ii] for ii in Op[i]) for jj in range(2)),\
                              sum(eigv[j][1][jj]*eigv[j][0][jj]*sum(bs1[jj][ii] for ii in Op[i]) for jj in range(2)),\
                              sum(eigv[j][1][jj]*eigv[j][1][jj]*sum(bs1[jj][ii] for ii in Op[i]) for jj in range(2))])
   return Op_val
          

def Construct_transition_elements(N_states,basis,eigv):
   tran_elm=[]
   tran_mat=[]
   size_elm=[]
   for j,bs1 in enumerate(basis):
      tran_elm.append([])
      tran_mat.append([])
      size_elm.append(len(bs1))
   for j,bs1 in enumerate(basis):
      for i in range(N_states):
         check=0
         if bs1[0][i]==0:
            bas1=copy.deepcopy(bs1[0]);check=1;n1=0
         elif (len(bs1)==2 and bs1[1][i]==0):
            bas1=copy.deepcopy(bs1[1]);check=1;n1=1
         else:
            tran_elm[j].append(0);tran_mat[j].append([])
         if check==1:
            for k,bs2 in enumerate(basis):
               chk2=0
               if bs2[0][i]==1:
                  bas2=copy.deepcopy(bs2[0]);chk2=1;n2=0
               elif (len(bs2)==2 and bs2[1][i]==1):
                  bas2=copy.deepcopy(bs2[1]);chk2=1;n2=1
               if chk2==1:
                  bass1=copy.deepcopy(bas1)
                  bass2=copy.deepcopy(bas2)
                  bass1.pop(i); bass2.pop(i)
                  if sum((array(bass1)-array(bass2))**2)==0:
                     tran_elm[j].append(k+1)
                     if len(bs1)==1 and len(bs2)==1:
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])])
                     elif len(bs1)==2:
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])*eigv[j][0][n1],(-1.0)**sum(bas1[:i])*eigv[j][1][n1]])
                     elif len(bs2)==2:
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])*eigv[k][0][n2],(-1.0)**sum(bas1[:i])*eigv[k][1][n2]])
                     else: print "Something is wrong!"; exit()
                     break
   return size_elm,tran_elm,tran_mat

def Construct_transition_elements_old(N_states,basis):
   tran_elm=[]
   tran_mat=[]
   size_elm=[]
   for j,bs1 in enumerate(basis):
      tran_elm.append([])
      tran_mat.append([])
      size_elm.append(len(bs1))
   for j,bs1 in enumerate(basis):
      for i in range(N_states):
         check=0
         if bs1[0][i]==0: 
            bas1=copy.deepcopy(bs1[0]);check=1;n1=0
         elif (len(bs1)==2 and bs1[1][i]==0):
            bas1=copy.deepcopy(bs1[1]);check=1;n1=1
         else: 
            tran_elm[j].append(0);tran_mat[j].append([])
         if check==1:
            for k,bs2 in enumerate(basis):
               chk2=0
               if bs2[0][i]==1:
                  bas2=copy.deepcopy(bs2[0]);chk2=1;n2=0
               elif (len(bs2)==2 and bs2[1][i]==1):
                  bas2=copy.deepcopy(bs2[1]);chk2=1;n2=1
               if chk2==1:
                  bass1=copy.deepcopy(bas1)
                  bass2=copy.deepcopy(bas2)
                  bass1.pop(i); bass2.pop(i)
                  if sum((array(bass1)-array(bass2))**2)==0:
                     tran_elm[j].append(k+1)
                     if len(bs1)==1 and len(bs2)==1:
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])])
                     elif n1==0 and n2==0:
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])/sqrt(2),(-1.0)**sum(bas1[:i])/sqrt(2)])
                     elif (n1==1 and n2==0) or (n1==0 and n2==1):
                        tran_mat[j].append([(-1.0)**sum(bas1[:i])/sqrt(2),(-1.0)**(sum(bas1[:i])+1)/sqrt(2)])
                     else: print "Something is wrong!"; exit()
                     break
#                     for ii in range(len(bs1)):
#                        for jj in range(len(bs2)):
#                           tran_mat[j].append([(-1.0)**sum(bas1[:i])*eigv[j][ii]*eigv[k][jj])
#               tran_elm[j].append([k+1,[eigv[(-1.0)**sum(bs1[:i])]]); check=1; break
#            if check==0: print "Could not find the matrix elements!"; exit()
#         elif len(bs1)==2 and bs1[1][i]==0:
   return size_elm,tran_elm,tran_mat


def Construct_matrix_elements(n_orb,n_spin,basis):
   N_states=n_orb*n_spin
   mat_elm=[]
   for i in range(N_states):
      mat_elm.append([])
      for j,bs1 in enumerate(basis):
         if bs1[i]==0:
            for k,bs2 in enumerate(basis):
               if bs2[i]==1:
                  bas1=copy.deepcopy(bs1); bas2=copy.deepcopy(bs2)
                  bas1.pop(i); bas2.pop(i)
                  if sum((array(bas1)-array(bas2))**2)==0:
                     mat_elm[i].append([[j,k],(-1)**sum(bs1[:i])]); break

   return mat_elm

def Print_cix(J,ed1,ed2):

   norb=5
   nspin=2
   N_states=norb*nspin
   sym=[0,0,0,0,1,1,1,1,1,1]
   Op=[[0,1,2,3],[4,5,6,7,8,9]]


   e_d=[0.0, 0.0, 0.0, 0.0, ed2-ed1, ed2-ed1, ed2-ed1, ed2-ed1, ed2-ed1, ed2-ed1]
   e_d2=[0.0, ed2-ed1]

   basis=Construct_basis(N_states)
   basis_pair=[]
   new_basis=[]
   for bs in basis:
      #idx=0
      #for bs2 in basis_pair:
      #   if bs==bs2[0]: new_basis.append([bs2[0],bs2[1]])
      #   elif bs!=bs2[1]: idx+=1
      #if idx==2: new_basis.append([bs])
      new_basis.append([bs])
   #print new_basis

   (E_atom,n_atom,Sz)=Construct_atomic_energy(J,e_d,norb,nspin,new_basis)
   #print Sz
   #print E_atom
   #print n_atom

   diag_E_atom,eigv=Off_J(new_basis,basis_pair,E_atom)
   #print 'diag_E_atom=',diag_E_atom
   #print 'eigv=',eigv

   Op_val=Construct_Operator(new_basis,eigv,Op)
   #print Op_val

#   size_elm,tran_elm,tran_mat=Construct_transition_elements_old(N_states,new_basis)
   size_elm,tran_elm,tran_mat=Construct_transition_elements(N_states,new_basis,eigv)
   #print tran_elm
   #print tran_mat

   fcix = open('impurity.cix', 'w')
   # ---------------------------------------------------------------------------------------
   # -------------- Below is printing for ctqmc  solver ------------------------------------
   # ---------------------------------------------------------------------------------------
   print >> fcix, '# CIX file for ctqmc! '
   print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# baths, dimension, symmetry'

   for ib in range(N_states):
      print >> fcix, ib, '  ', 1, sym[ib], sym[ib] 

   print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
   for E in e_d2:  print >> fcix, E,
   print >> fcix
   print >> fcix, '#   N   K   Sz size'


   for i in range(len(new_basis)):
      print >> fcix, "%3d  %2d %2d %4.1f %2d " % (i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""   
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else: 
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

   print >> fcix, 'HB1'
   print >> fcix, '# number of operators needed'
   print >> fcix, 1
   print >> fcix, '# Occupancy'
   for ib in range(len(new_basis)):
      for jb in range(len(Op_val)):
         print >> fcix, "%3d %3d  %2d " % (ib+1, int(sqrt(len(Op_val[jb][ib]))), int(sqrt(len(Op_val[jb][ib])))), 
         for kb in range(len(Op_val[jb][ib])): print >> fcix, Op_val[jb][ib][kb],
         print >> fcix, ''
   print >> fcix, '# Data for HB1'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# ind   N   K   Jz size'
   for i in range(len(new_basis)):
      print >> fcix, "%3d %3d  %2d %2d %4.1f %2d " % (i+1, i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else:
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

def Print_cix_3orb(J,ed1,ed2,ed3):

   norb=5
   nspin=2
   N_states=norb*nspin
   sym=[0,0,1,1,1,1,2,2,2,2]
   Op=[[0,1],[2,3,4,5],[6,7,8,9]]


   e_d=[0.0, 0.0, ed2-ed1, ed2-ed1, ed2-ed1, ed2-ed1, ed3-ed1, ed3-ed1, ed3-ed1, ed3-ed1]
   e_d2=[0.0, ed2-ed1, ed3-ed1]

   basis=Construct_basis(N_states)
   basis_pair=[]
   new_basis=[]
   for bs in basis:
      #idx=0
      #for bs2 in basis_pair:
      #   if bs==bs2[0]: new_basis.append([bs2[0],bs2[1]])
      #   elif bs!=bs2[1]: idx+=1
      #if idx==2: new_basis.append([bs])
      new_basis.append([bs])
   #print new_basis

   (E_atom,n_atom,Sz)=Construct_atomic_energy(J,e_d,norb,nspin,new_basis)
   #print Sz
   #print E_atom
   #print n_atom

   diag_E_atom,eigv=Off_J(new_basis,basis_pair,E_atom)
   #print 'diag_E_atom=',diag_E_atom
   #print 'eigv=',eigv

   Op_val=Construct_Operator(new_basis,eigv,Op)
   #print Op_val

#   size_elm,tran_elm,tran_mat=Construct_transition_elements_old(N_states,new_basis)
   size_elm,tran_elm,tran_mat=Construct_transition_elements(N_states,new_basis,eigv)
   #print tran_elm
   #print tran_mat

   fcix = open('impurity.cix', 'w')
   # ---------------------------------------------------------------------------------------
   # -------------- Below is printing for ctqmc  solver ------------------------------------
   # ---------------------------------------------------------------------------------------
   print >> fcix, '# CIX file for ctqmc! '
   print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# baths, dimension, symmetry'

   for ib in range(N_states):
      print >> fcix, ib, '  ', 1, sym[ib], sym[ib] 

   print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
   for E in e_d2:  print >> fcix, E,
   print >> fcix
   print >> fcix, '#   N   K   Sz size'


   for i in range(len(new_basis)):
      print >> fcix, "%3d  %2d %2d %4.1f %2d " % (i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""   
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else: 
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

   print >> fcix, 'HB1'
   print >> fcix, '# number of operators needed'
   print >> fcix, 1
   print >> fcix, '# Occupancy'
   for ib in range(len(new_basis)):
      for jb in range(len(Op_val)):
         print >> fcix, "%3d %3d  %2d " % (ib+1, int(sqrt(len(Op_val[jb][ib]))), int(sqrt(len(Op_val[jb][ib])))), 
         for kb in range(len(Op_val[jb][ib])): print >> fcix, Op_val[jb][ib][kb],
         print >> fcix, ''
   print >> fcix, '# Data for HB1'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# ind   N   K   Jz size'
   for i in range(len(new_basis)):
      print >> fcix, "%3d %3d  %2d %2d %4.1f %2d " % (i+1, i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else:
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""


def Print_cix_5orb(Ueff,ed):#1,ed2,ed3,ed4,ed5):

   norb=5
   nspin=2
   N_states=norb*nspin
   sym=[0,0,1,1,2,2,3,3,4,4]
   Op=[[0,1],[2,3],[4,5],[6,7],[8,9]]

   ed1=ed[0];ed2=ed[1];ed3=ed[2];ed4=ed[3];ed5=ed[4]
   e_d=[0.0, 0.0, ed2-ed1, ed2-ed1, ed3-ed1, ed3-ed1, ed4-ed1, ed4-ed1, ed5-ed1, ed5-ed1]
   e_d2=[0.0, ed2-ed1, ed3-ed1, ed4-ed1, ed5-ed1]

   basis=Construct_basis(N_states)
   basis_pair=[]
   new_basis=[]
   for bs in basis:
      #idx=0
      #for bs2 in basis_pair:
      #   if bs==bs2[0]: new_basis.append([bs2[0],bs2[1]])
      #   elif bs!=bs2[1]: idx+=1
      #if idx==2: new_basis.append([bs])
      new_basis.append([bs])
   #print new_basis

   (E_atom,n_atom,Sz)=Construct_atomic_energy2(Ueff,e_d,norb,nspin,new_basis)
   #print Sz
   #print E_atom
   #print n_atom

   diag_E_atom,eigv=Off_J(new_basis,basis_pair,E_atom)
   #print 'diag_E_atom=',diag_E_atom
   #print 'eigv=',eigv

   Op_val=Construct_Operator(new_basis,eigv,Op)
   #print Op_val

#   size_elm,tran_elm,tran_mat=Construct_transition_elements_old(N_states,new_basis)
   size_elm,tran_elm,tran_mat=Construct_transition_elements(N_states,new_basis,eigv)
   #print tran_elm
   #print tran_mat

   fcix = open('impurity.cix', 'w')
   # ---------------------------------------------------------------------------------------
   # -------------- Below is printing for ctqmc  solver ------------------------------------
   # ---------------------------------------------------------------------------------------
   print >> fcix, '# CIX file for ctqmc! '
   print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# baths, dimension, symmetry'

   for ib in range(N_states):
      print >> fcix, ib, '  ', 1, sym[ib], sym[ib] 

   print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
   for E in e_d2:  print >> fcix, E,
   print >> fcix
   print >> fcix, '#   N   K   Sz size'


   for i in range(len(new_basis)):
      print >> fcix, "%3d  %2d %2d %4.1f %2d " % (i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""   
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else: 
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

   print >> fcix, 'HB1'
   print >> fcix, '# number of operators needed'
   print >> fcix, 1
   print >> fcix, '# Occupancy'
   for ib in range(len(new_basis)):
      for jb in range(len(Op_val)):
         print >> fcix, "%3d %3d  %2d " % (ib+1, int(sqrt(len(Op_val[jb][ib]))), int(sqrt(len(Op_val[jb][ib])))), 
         for kb in range(len(Op_val[jb][ib])): print >> fcix, Op_val[jb][ib][kb],
         print >> fcix, ''
   print >> fcix, '# Data for HB1'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# ind   N   K   Jz size'
   for i in range(len(new_basis)):
      print >> fcix, "%3d %3d  %2d %2d %4.1f %2d " % (i+1, i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else:
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""


def Print_cix_sp(J,ed1,ed2,ed3,esplit):

   norb=5
   nspin=2
   N_states=norb*nspin
   sym=[0,1,2,3,2,3,4,5,4,5]
   Op=[[0],[1],[2,4],[3,5],[6,8],[7,9]]


   e_d=[0.0, esplit, ed2-ed1, ed2-ed1+esplit, ed2-ed1, ed2-ed1+esplit, ed3-ed1, ed3-ed1+esplit, ed3-ed1, ed3-ed1+esplit]
   e_d2=[0.0, esplit, ed2-ed1, ed2-ed1+esplit, ed3-ed1, ed3-ed1+esplit]

   basis=Construct_basis(N_states)
   basis_pair=[]
   new_basis=[]
   for bs in basis:
      #idx=0
      #for bs2 in basis_pair:
      #   if bs==bs2[0]: new_basis.append([bs2[0],bs2[1]])
      #   elif bs!=bs2[1]: idx+=1
      #if idx==2: new_basis.append([bs])
      new_basis.append([bs])
   #print new_basis

   (E_atom,n_atom,Sz)=Construct_atomic_energy(J,e_d,norb,nspin,new_basis)
   #print Sz
   #print E_atom
   #print n_atom

   diag_E_atom,eigv=Off_J(new_basis,basis_pair,E_atom)
   #print 'diag_E_atom=',diag_E_atom
   #print 'eigv=',eigv

   Op_val=Construct_Operator(new_basis,eigv,Op)
   #print Op_val

#   size_elm,tran_elm,tran_mat=Construct_transition_elements_old(N_states,new_basis)
   size_elm,tran_elm,tran_mat=Construct_transition_elements(N_states,new_basis,eigv)
   #print tran_elm
   #print tran_mat

   fcix = open('impurity.cix', 'w')
   # ---------------------------------------------------------------------------------------
   # -------------- Below is printing for ctqmc  solver ------------------------------------
   # ---------------------------------------------------------------------------------------
   print >> fcix, '# CIX file for ctqmc! '
   print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# baths, dimension, symmetry'

   for ib in range(N_states):
      print >> fcix, ib, '  ', 1, sym[ib], sym[ib] 

   print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
   for E in e_d2:  print >> fcix, E,
   print >> fcix
   print >> fcix, '#   N   K   Sz size'


   for i in range(len(new_basis)):
      print >> fcix, "%3d  %2d %2d %4.1f %2d " % (i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""   
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else: 
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

   print >> fcix, 'HB1'
   print >> fcix, '# number of operators needed'
   print >> fcix, 1
   print >> fcix, '# Occupancy'
   for ib in range(len(new_basis)):
      for jb in range(len(Op_val)):
         print >> fcix, "%3d %3d  %2d " % (ib+1, int(sqrt(len(Op_val[jb][ib]))), int(sqrt(len(Op_val[jb][ib])))), 
         for kb in range(len(Op_val[jb][ib])): print >> fcix, Op_val[jb][ib][kb],
         print >> fcix, ''
   print >> fcix, '# Data for HB1'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# ind   N   K   Jz size'
   for i in range(len(new_basis)):
      print >> fcix, "%3d %3d  %2d %2d %4.1f %2d " % (i+1, i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else:
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""


def Print_cix_sp_5orb(J,ed1,ed2,ed3,ed4,ed5,esplit):

   norb=5
   nspin=2
   N_states=norb*nspin
   sym=[0,1,2,3,4,5,6,7,8,9]
   Op=[[0],[1],[2],[3],[4],[5],[6],[7],[8],[9]]


   e_d=[0.0, esplit, ed2-ed1, ed2-ed1+esplit, ed3-ed1, ed3-ed1+esplit, ed4-ed1, ed4-ed1+esplit, ed5-ed1, ed5-ed1+esplit]
   e_d2=[0.0, esplit, ed2-ed1, ed2-ed1+esplit, ed3-ed1, ed3-ed1+esplit, ed4-ed1, ed4-ed1+esplit, ed5-ed1, ed5-ed1+esplit]
   #e_d2=[0.0, esplit, ed2-ed1, ed2-ed1+esplit, ed3-ed1, ed3-ed1+esplit]

   basis=Construct_basis(N_states)
   basis_pair=[]
   new_basis=[]
   for bs in basis:
      #idx=0
      #for bs2 in basis_pair:
      #   if bs==bs2[0]: new_basis.append([bs2[0],bs2[1]])
      #   elif bs!=bs2[1]: idx+=1
      #if idx==2: new_basis.append([bs])
      new_basis.append([bs])
   #print new_basis

   (E_atom,n_atom,Sz)=Construct_atomic_energy(J,e_d,norb,nspin,new_basis)
   #print Sz
   #print E_atom
   #print n_atom

   diag_E_atom,eigv=Off_J(new_basis,basis_pair,E_atom)
   #print 'diag_E_atom=',diag_E_atom
   #print 'eigv=',eigv

   Op_val=Construct_Operator(new_basis,eigv,Op)
   #print Op_val

#   size_elm,tran_elm,tran_mat=Construct_transition_elements_old(N_states,new_basis)
   size_elm,tran_elm,tran_mat=Construct_transition_elements(N_states,new_basis,eigv)
   #print tran_elm
   #print tran_mat

   fcix = open('impurity.cix', 'w')
   # ---------------------------------------------------------------------------------------
   # -------------- Below is printing for ctqmc  solver ------------------------------------
   # ---------------------------------------------------------------------------------------
   print >> fcix, '# CIX file for ctqmc! '
   print >> fcix, '# cluster_size, number of states, number of baths, maximum_matrix_size'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# baths, dimension, symmetry'

   for ib in range(N_states):
      print >> fcix, ib, '  ', 1, sym[ib], sym[ib] 

   print >> fcix, '# cluster energies for non-equivalent baths, eps[k]'
   for E in e_d2:  print >> fcix, E,
   print >> fcix
   print >> fcix, '#   N   K   Sz size'


   for i in range(len(new_basis)):
      print >> fcix, "%3d  %2d %2d %4.1f %2d " % (i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""   
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else: 
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""

   print >> fcix, 'HB1'
   print >> fcix, '# number of operators needed'
   print >> fcix, 1
   print >> fcix, '# Occupancy'
   for ib in range(len(new_basis)):
      for jb in range(len(Op_val)):
         print >> fcix, "%3d %3d  %2d " % (ib+1, int(sqrt(len(Op_val[jb][ib]))), int(sqrt(len(Op_val[jb][ib])))), 
         for kb in range(len(Op_val[jb][ib])): print >> fcix, Op_val[jb][ib][kb],
         print >> fcix, ''
   print >> fcix, '# Data for HB1'
   print >> fcix, 1, len(new_basis), N_states, max(size_elm)
   print >> fcix, '# ind   N   K   Jz size'
   for i in range(len(new_basis)):
      print >> fcix, "%3d %3d  %2d %2d %4.1f %2d " % (i+1, i+1, n_atom[i][0], 0, Sz[i][0], size_elm[i]),#len(Enes[ii])),
      for ib in range(N_states):
         print >> fcix, "%3d" % (tran_elm[i][ib]),
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, diag_E_atom[i][j]," ",
      print >> fcix, "  ",
      for j in range(size_elm[i]): print >> fcix, abs(Sz[j][0])," ",
      print >> fcix, ""
   print >> fcix, '# matrix elements'
   for i,elm in enumerate(tran_elm):
      for ib in range(N_states):
         print >> fcix, "%3d %3d "  % (i+1, tran_elm[i][ib]),
         if tran_elm[i][ib]==0:
            print >> fcix, "%2d %2d" % (0, 0)
         else:
            print >> fcix, "%2d %2d " % (size_elm[i], size_elm[tran_elm[i][ib]-1]),
            for k in range(size_elm[i]*size_elm[tran_elm[i][ib]-1]): print >> fcix, tran_mat[i][ib][k],
            print >> fcix, ""
