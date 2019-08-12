#!/usr/bin/env python
import sys, subprocess, os
import numpy as np
from scipy import *
import copy, Fileio, re
from scipy import interpolate
import shutil
from shutil import copyfile

#######inputs######################################################################
#dmft_bin
path_bin="/home/uthpala/Documents/Research/projects/DMFTwDFT/bin/"

#mpirun
if os.path.exists("para_com.dat"):
    fipa=open('para_com.dat','r')
    para_com=str(fipa.readline())[:-1]
    fipa.close()
else:
    para_com=""

#DOS parameters
emin=-7.0
emax=3.0
rom=3000 #omega points
broaden=0.03

###################################################################################

print('\n###############################')
print('# DMFT post-processing scheme #')
print('###############################\n')

#creating directory for dos
if os.path.exists("dos"):
	shutil.rmtree("dos")
	os.makedirs("dos")
else:	
	os.makedirs("dos")

#copying the last few self-energies from the DMFT run in the directory above	
os.popen("cp -r ../sig.inp.0.* ./")

#averaging sef energies
print('Averaging self-energies...')
cmd = "sigaver.py"
out, err = subprocess.Popen(cmd, shell=True).communicate()
print('Complete.\n')

#copy maxent_params.dat from source
src=path_bin+ os.sep+"maxent_params.dat"
copyfile(src,"./")

#Analytic continuation
print('Analytic Continuation...\n')
cmd = "maxent_run.py sig.inpx"
subprocess.Popen(cmd, shell=True).communicate()
print('Complete.\n')

#copying files from DMFT directory
cmd = "Copy_input.py -dos ../"
subprocess.Popen(cmd, shell=True).communicate()

#interpolating points on real axis
headerline=2
om,Sig=Fileio.Read_complex_multilines('Sig.out',headerline)
s_oo = None
Vdc = None
fi=open('Sig.out','r')
for i in range(headerline):
  line=fi.readline()
  m=re.search('#(.*)',line)
  exec(m.group(1).strip())
s_oo_Vdc=array(s_oo)-array(Vdc)

ommesh=linspace(emin,emax,rom)
Sig_tot=zeros((len(Sig),rom),dtype=complex)
for i in range(len(Sig)):
  SigSpline = interpolate.splrep(om, Sig[i].real, k=1, s=0)
  Sig_tot[i,:] += interpolate.splev(ommesh, SigSpline)
  SigSpline = interpolate.splrep(om, Sig[i].imag, k=1, s=0)
  Sig_tot[i,:] += 1j*interpolate.splev(ommesh, SigSpline)

header1='# nom,ncor_orb= '+str(len(ommesh))+' '+str(len(Sig_tot))
#header2='# T= %18.15f'%(1.0/pC['beta'][0])#+str(self.T)
header2='# T= %18.15f'%(broaden)#+str(self.T)
header3='# s_oo-Vdc= '
for i in range(len(s_oo_Vdc)):
  header3+='%18.15f '%(s_oo_Vdc[i])
header4='# s_oo= '+str(s_oo)
header5='# Vdc= '+str(Vdc)
Fileio.Print_complex_multilines(Sig_tot,ommesh,'sig.inp_real',[header1,header2,header3,header4,header5])

#running dmft_dos.x
cmd = para_com + "dmft_dos.x"
subprocess.Popen(cmd, shell=True).communicate()
