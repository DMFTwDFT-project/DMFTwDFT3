#!/usr/bin/env python2
from scipy import *
import sys, copy, re

def Read_complex_multilines(D_name,skipline=0):
   """ This function reads Sig.out file """
   print("Reading a file ",D_name)
   fi=open(D_name,'r')
   for i in range(skipline):
      fi.readline()
   lines=fi.readlines()
   fi.close()
   m=len(lines[0].split())
   Nom=len(lines)
   om_data=zeros(Nom,dtype=float)
   Data=zeros(((m-1)/2,Nom), dtype=complex)
   for iom,line in enumerate(lines):
      l=line.split()
      om_data[iom]=float(l[0])
      for ib in range((m-1)/2):
         Data[ib,iom]=complex(float(l[1+ib*2]),float(l[2+ib*2]))
   return om_data, Data

def Read_float_multilines(D_name):
   """ This function reads Sig.out file """
   print("Reading a file ",D_name)
   fi=open(D_name,'r')
   lines=fi.readlines()
   fi.close()
   m=len(lines[0].split())
   Nom=len(lines)
   #om_data=zeros(Nom,dtype=float)
   Data=zeros((m,Nom))
   for iom,line in enumerate(lines):
      l=line.split()
      for ib in range(m):
         Data[ib,iom]=float(l[ib])
   return Data

def Read_float(D_name):
   """ This function reads Sig.out file """
   print("Reading a file ",D_name)
   fi=open(D_name,'r')
   lines=fi.readlines()
   fi.close()
   Data=array(list(map(float,lines[0].split())))
   return Data

def Read_complex_Data(D_name):
   """ This function reads Sig.out file """
   print("Reading a file ",D_name)
   fi=open(D_name,'r')
   line=fi.readline()
   val=re.search(r'TrSigmaG=(\-?\d+\.?\d*)',line)
   TrSigmaG=float(val.group(1))
   val=re.search(r'mu=(\-?\d+\.?\d*)',line)
   mu=float(val.group(1))
   val=re.search(r'Ekin=(\-?\d+\.?\d*)',line)
   Ekin=float(val.group(1))
   val=re.search(r'Epot=(\-?\d+\.?\d*)',line)
   Epot=float(val.group(1))
   val=re.search(r'nf=(\d+\.?\d*)',line)
   nf_q=float(val.group(1))
   mom=eval(line.split()[-1][4:])
   lines=fi.readlines()
   Nom=len(lines)
   m=len(lines[0].split())
   om_data=zeros(Nom,dtype=float)
   Data=zeros(((m-1)/2,Nom), dtype=complex)
   iom=0
   for line in lines:
      l=line.split()
      om_data[iom]=float(l[0])
      for ib in range((m-1)/2):
         Data[ib,iom]=complex(float(l[1+ib*2]),float(l[2+ib*2]))
      iom=iom+1
   if iom!=Nom:
      print("Something is wrong!")
      sys.exit(1)
   if len(mom)!=(m-1)/2:
      print("Something is wrong!")
      sys.exit(1)

   return (om_data,Data,TrSigmaG,Epot,nf_q,mom,Ekin,mu)

def Print_complex(data,mesh,filename):
   n1=len(mesh)
   fi=open(filename,'w')
   for i in range(n1):
      f.write("%.14f " %(mesh[i]))
      f.write("%.14f %.14f " %(data[i].real, data[i].imag))
   fi.close()

def Print_complex_multilines(data,mesh,filename,headers=[]):
   n0=len(data)
   n1=len(mesh)
   fi=open(filename,'w')
   for header in headers:
      fi.write("$s" % header)
   for i in range(n1):
      for j in range(n0):
         if j==0: fi.write("%20.15f " %(mesh[i]))
         fi.write("%20.15f %20.15f " %(data[j,i].real, data[j,i].imag))
      fi.write("")
   fi.close()

def Print_float(data,filename):
   n0=len(data)
   fi=open(filename,'w')
   for j in range(n0):
      fi.write("%.14f " %(data[j]))
   fi.write("")
   fi.close()

def Read_float(filename):
   fi=open(filename,'r')
   Data=array(list(map(float,fi.readline().split())))
   return Data

def Create_dmft_params(p,pC,N_atoms,atm_idx,sym_idx):
   f = open('dmft_params.dat', 'w')
   f.write("# Number of k-points in Wannier basis=\n")
   f.write("%d %d %d\n " % (p['q'][0],p['q'][1],p['q'][2]))
   f.write("# Total number of electrons=\n")
   f.write("%d\n" % p['n_tot'])
#   print >> f, "# Temperature [eV]="
#   print >> f, 1.0/pC['beta'][0]
   f.write("# Number of om points for k-sum\n")
   f.write("%d\n" % p['noms'])
   f.write("# Number of iterations for mu\n")
   f.write("%d\n" % p['mu_iter'])
   f.write("# Number of total spin\n")
   f.write("%d\n" % p['nspin'])
   f.write("# Number of total correlated atoms\n")
   f.write("%d\n" % N_atoms)
   f.write("# Number of correlated orbitals per atom\n")
   f.write("%d\n" % len(sym_idx[atm_idx[0]]))
   f.write("# Orbital index for the self-energy at each atom\n")
   for i in range(N_atoms):
      for j in range(len(sym_idx[atm_idx[i]])):
         f.write("%d " % sym_idx[atm_idx[i]][j])
      f.write(" ")
   f.close()

def Create_INPUT(p,pC,TB,T_high,noms_high,LFORCE=".FALSE."):
   atm_idx=[]
   idx=1
   for ats in p['cor_at']:
      for at in ats:
         atm_idx.append(idx)
      idx+=1
   f = open('VASP.input', 'w')
   f.write("%s" % LFORCE)    
   f.write("%s" % TB.LHF)
   f.write("%d" % p['n_tot'])
   f.write("%d" % p['nspin'])
   f.write("%d" % p['nfine'])
   f.write("%s" % TB.ncor_orb)
   f.write("%s" % TB.max_cor_orb)
   if TB.LHF==".TRUE.":
      f.write('1')
   else:
      f.write("%d" % p['noms'])
   if TB.LHF==".TRUE.":
      f.write('1')
   else:
      f.write("%f" % noms_high) 
   if TB.LHF==".TRUE.":
      f.write('1')
   else:
      f.write("%d " % (p['noms']+p['nomlog']))
   f.write("%f" % (1.0/pC['beta'][0]))
   f.write("%f" % T_high) 
   for i in range(len(atm_idx)):
      f.write("%s " % atm_idx[i])
   f.write('')
   for i in range(len(atm_idx)):
      f.write("%f" % p['U'][atm_idx[i]-1])
   f.write('')
   for i in range(len(atm_idx)):
      f.write("%f" % p['J'][atm_idx[i]-1])
   f.write('')
   for i in range(len(atm_idx)):
      f.write("%f" % p['alpha'][atm_idx[i]-1])
   f.write('')
   f.close()
   if TB.LHF==".FALSE.":
      f = open('ksum.input', 'w')
      f.write("%d %d %d " % (p['q'][0],p['q'][1],p['q'][2]))
      f.write("%d %d" % (['noms'], (p['noms']+p['nomlog'])))
      f.write("%d" % p['nspin']) 
      f.write("%d" % TB.ncor_orb) 
      f.write("%d" % TB.max_cor_orb)
      for i in range(len(atm_idx)): 
         f.write("%s" % atm_idx[i])
      f.write('')
      f.write("%f" % 1.0/pC['beta'][0])
#      print >> f, p['Nd_f'] 
      f.write("%d" % p['n_tot'])
      f.write("%d" % p['mu_iter'])
      f.write("%f" % p['mix_sig'])
      for i in range(len(atm_idx)):
         f.write("%f" % p['U'][atm_idx[i]-1])
      f.write('') 
      for i in range(len(atm_idx)):
         f.write("%f" % p['alpha'][atm_idx[i]-1])
      f.write('') 
      for i in range(len(atm_idx)):
         f.write("%f" % p['J'][atm_idx[i]-1])
      f.write('') 
      f.close()
   else:
      f = open('ksum.input', 'w')
      f.write("%d" % p['nspin'])
      f.write("%d" % TB.ncor_orb)
      f.write("%d" % TB.max_cor_orb)
      f.write("%d" % p['n_tot'])
      f.write("%d" % p['mu_iter'])
      for i in range(len(atm_idx)):
         f.write("%s" % atm_idx[i])
      f.write('')
      for i in range(len(atm_idx)):
         f.write("%f" % p['U'][atm_idx[i]-1])
      f.write('')
      for i in range(len(atm_idx)):
         f.write("%f" % p['alpha'][atm_idx[i]-1])
      f.write('')
      for i in range(len(atm_idx)):
         f.write("%f" % p['J'][atm_idx[i]-1])
      f.write('')
      f.close()
