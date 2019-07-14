#!/usr/bin/env python2
import sys, subprocess, os
import numpy as np
import shutil
from shutil import copyfile
import VASP
import Struct

###############################inputs######################################################################
#number of bands for DFT calculation
NUM_BANDS = 72

#vasp executable
vasp_exec = "vasp_std"

#dictionary of values for generating wannier90.win
gen_win_dic = {
	"atomnames": ['V','O'],       # The name of atoms"],
	"orbs"     : ['d','p'],       # The name  of orbitals
	"L_rot"     : [1,0],       # Whether rotate local axis or not
	"ewin":      [-7,7],           # Energy Window
    "cor_at":    [['V1']],      # Correlated atoms, put degenerate atoms in the same list
    "cor_orb":   [[['d_z2','d_x2y2'],['d_xz','d_yz','d_xy']]], # DMFT orbitals, other orbitals are treated by HF"],
}


#mpirun
if os.path.exists("para_com.dat"):
    fipa=open('para_com.dat','r')
    para_com=str(fipa.readline())[:-1]
    fipa.close()
else:
    para_com=""

###################################################################################

print('\n#######################')
print('# DMFTwDFT initialization #')
print('##########################\n')

############initialization############################################

#generating wannier90.win
TB=Struct.TBstructure('POSCAR',gen_win_dic['atomnames'],gen_win_dic['orbs'])
TB.Compute_cor_idx(gen_win_dic['cor_at'],gen_win_dic['cor_orb'])
print(TB.TB_orbs)
DFT=VASP.VASP_class()
DFT.NBANDS=NUM_BANDS
DFT.Create_win(TB,gen_win_dic['atomnames'],gen_win_dic['orbs'],gen_win_dic['L_rot'],DFT.NBANDS,DFT.EFERMI+gen_win_dic['ewin'][0],DFT.EFERMI+gen_win_dic['ewin'][1])

#initial DFT run
print('Running initial DFT...')
cmd = para_com+" "+vasp_exec
out, err = subprocess.Popen(cmd, shell=True).communicate()
print('Initial DFT calculation complete.\n')

#running wannier90.x to generate .chk
print('Running wannier90...')
cmd = "wannier90.x wannier90"
out, err = subprocess.Popen(cmd, shell=True).communicate()
print('wannier90 calculation complete.\n')

#generate sig.inp 
cmd = "sigzero.py"
out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
print('Generated sig.inp.\n')

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
print('DMFT initialization complete.\n')

