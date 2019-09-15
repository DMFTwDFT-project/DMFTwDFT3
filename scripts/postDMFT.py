#!/usr/bin/env python
import sys, subprocess, os
import numpy as np
from scipy import *
import copy, Fileio, re
from scipy import interpolate
import shutil
import glob
from INPUT import *

#######inputs######################################################################

#DOS parameters
emin=-7.0
emax=3.0
rom=3000 #omega points
broaden=0.03

#How many last self energies?
siglistindx=5

###################################################################################

#mpirun
if os.path.exists("para_com.dat"):
    fipa=open('para_com.dat','r')
    para_com=str(fipa.readline())[:-1]
    fipa.close()
else:
    para_com=""


print('\n')
print('-------------------------------')
print('| DMFT post-processing scheme |')
print('-------------------------------\n')

#dmft_bin
path_bin=p['path_bin']

#creating directory for dos
if os.path.exists("dos"):
	shutil.rmtree("dos")
	os.makedirs("dos")
else:	
	os.makedirs("dos")

#copying the last few self-energies from the DMFT run in the directory above
siglist = sorted(glob.glob("sig.inp.*"),key=os.path.getmtime)[-siglistindx:]
for file in siglist:
	shutil.copy(file,'./dos')

#averaging sef energies
print('Averaging self-energies from: ')
print(siglist)
cmd = "cd dos && sigaver.py sig.inp.*"
out, err = subprocess.Popen(cmd, shell=True,stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
if err:
  print(err)
  print('Averaging self-energies Failed!\n')
  sys.exit()
else:
  print('Self-energies averaged.\n')  
  

#copy maxent_params.dat from source if not in DMFT directory
if os.path.exists("maxent_params.dat"):
  shutil.copyfile('maxent_params.dat','./dos/maxent_params.dat')
else:
  src=path_bin+ os.sep+"maxent_params.dat"
  shutil.copyfile(src,"./dos/maxent_params.dat")

#Analytic continuation
print('Running analytic continuation...')
cmd = "cd dos && maxent_run.py sig.inpx"+" > ac.out 2> ac.error"
out, err = subprocess.Popen(cmd, shell=True,stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
if os.path.exists("./dos/Sig.out"):
  print('Analytic continuation complete.\n')  
else:
  print('Analytic continuation failed! Check ac.error for details.\n') 
  sys.exit()
  

#copying files from DMFT directory
cmd = "cd dos && Copy_input.py ../ -dos"
out, err = subprocess.Popen(cmd, shell=True).communicate()
if err:
  print('File copy failed!\n')
  print(err)
  sys.exit()
else:
  print(out)

#interpolating points on real axis
print('\nInterpolating points on real axis...') 
headerline=2
om,Sig=Fileio.Read_complex_multilines('./dos/Sig.out',headerline)
s_oo = None
Vdc = None
fi=open('./dos/Sig.out','r')
for i in range(headerline):
  line=fi.readline()
  m=re.search('#(.*)',line)
  exec(m.group(1).strip())
#s_oo_Vdc=array(s_oo)-array(Vdc)
s_oo_Vdc=array((np.array(s_oo)).astype(np.float))-array((np.array(Vdc)).astype(np.float))

ommesh=linspace(emin,emax,rom)
Sig_tot=zeros((len(Sig),rom),dtype=complex)
for i in range(len(Sig)):
  SigSpline = interpolate.splrep(om, Sig[i].real, k=1, s=0)
  Sig_tot[i,:] += interpolate.splev(ommesh, SigSpline)
  SigSpline = interpolate.splrep(om, Sig[i].imag, k=1, s=0)
  Sig_tot[i,:] += 1j*interpolate.splev(ommesh, SigSpline)

fo=open('SigMdc.out','w')
for i in range(5):
   if i==0 or i==3:
      fo.write(%18.15f ' %(s_oo_Vdc[0]))
   else:
      fo.write(%18.15f ' %(s_oo_Vdc[1]))
	
header1='# nom,ncor_orb= '+str(len(ommesh))+' '+str(len(Sig_tot))
#header2='# T= %18.15f'%(1.0/pC['beta'][0])#+str(self.T)
header2='# T= %18.15f'%(broaden)#+str(self.T)
header3='# s_oo-Vdc= '
for i in range(len(s_oo_Vdc)):
  header3+='%18.15f '%(s_oo_Vdc[i])
header4='# s_oo= '+str(s_oo)
header5='# Vdc= '+str(Vdc)
Fileio.Print_complex_multilines(Sig_tot,ommesh,'./dos/sig.inp_real',[header1,header2,header3,header4,header5])
print('Interpolation complete.\n')

#running dmft_dos.x
print("Calculating DMFT DOS...")
cmd ="cd dos && "+para_com+" "+"dmft_dos.x"
out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
if err:
  print(err)
  print('DMFT DOS calculation failed!\n')
  sys.exit()
else:
  print('DMFT DOS calculation complete.\n')  
print("DMFT post-processing complete")

