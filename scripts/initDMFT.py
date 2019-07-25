#!/usr/bin/env python2
import sys, subprocess, os
import numpy as np
import shutil
from shutil import copyfile
import VASP
import Struct
from INPUT import *


#vasp executable
vasp_exec = "vasp_std"

#mpirun
if os.path.exists("para_com.dat"):
    fipa=open('para_com.dat','r')
    para_com=str(fipa.readline())[:-1]
    fipa.close()
else:
    para_com=""

print('\n')
print('---------------------------')
print('| DMFTwDFT initialization |')
print('---------------------------\n')

############initialization############################################

#generating wannier90.win
TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
TB.Compute_cor_idx(p['cor_at'],p['cor_orb'])
print(TB.TB_orbs)
DFT=VASP.VASP_class()
DFT.NBANDS=pV['NBANDS='][0]
DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1])

#initial DFT run
print('\nRunning initial DFT...')
cmd = para_com+" "+vasp_exec #+ " > dft.out 2> dft.error"
out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
if err:
	print('DFT calculation failed! Check dft.error for details.\n')
	f = open('dft.error','w')
	f.write(err)
	f.close()
	sys.exit()
else:
	print('DFT calculation complete.\n')	
	f = open('dft.out','w')
	f.write(out)
	f.close()
	


#running wannier90.x to generate .chk
print('Running wannier90...')
cmd = para_com+" "+"wannier90.x wannier90"
out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
if err:
	print('wannier90 calculation failed!\n')
	print(err)
	sys.exit()
else:
	print('wannier90 calculation complete.\n')	
	print(out)


#generate sig.inp 
cmd = "sigzero.py"
out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
if err:
	print(err)
	sys.exit()
else:
	print('Initial self-energy file generated.\n')
	
#creating directory for DMFT
if os.path.exists("DMFT"):
	shutil.rmtree("DMFT")
	os.makedirs("DMFT")
else:	
	os.makedirs("DMFT")

#copy INPUT.py to DMFT directory
copyfile("INPUT.py","./DMFT/INPUT.py")

#copying files into DMFT directory
cmd = "cd ./DMFT && Copy_input.py ../"
out, err = subprocess.Popen(cmd, shell=True).communicate()
if err:
	print('File copy failed!\n')
	print(err)
	sys.exit()
else:
	print(out)
print('\nDMFT initialization complete. Ready to run RUNDMFT.py.\n')

