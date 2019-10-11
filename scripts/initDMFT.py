#!/usr/bin/env python3
import sys, subprocess, os
import numpy as np
import shutil
from shutil import copyfile
import VASP3
import Struct
from INPUT import *
import argparse
from argparse import RawTextHelpFormatter
import pymatgen.io.vasp as vasp


class Initialize():
	"""DMFTwDFT Initialization

	This class contains methods to run the initial DFT and wannier90 calculation
	to generate inputs for the DMFT runs.

	For ionic convergence create a directory DFT_relax in the root directory and run a convergence 
	calculation. Do this calculation within the submission script. 

	E.g.:

	cd $PBS_O_WORKDIR/DFT
	mpirun -np 40 vasp_std
	cp CONTCAR ../POSCAR
	cd ..
	python initDMFT.py
	cd DMFT
	RUNDMFT.py 

	"""

	def __init__(self,args):
		"""
		Contains common functions for all methods.
		"""
		if os.path.exists("para_com.dat"):
			fipa=open('para_com.dat','r')
			self.para_com=str(fipa.readline())[:-1]
			fipa.close()
		else:
			self.para_com=""

		#VASP	
		if args.dft == 'vasp':	
		
			#vasp executable
			self.vasp_exec = args.dftexec

			#import the VASP class
			self.DFT = VASP3.VASP_class()

			#vasp running directory (current directory)
			self.dir = os.getcwd()
					
		self.gen_win()	
		self.gen_sig()
		if args.relax:
			self.vasp_convergence()
		self.vasp_run(self.dir)
		self.update_win()
		self.run_wan90()
		self.copy_files()

	def gen_win(self):
		"""
		This method generates wannier90.win for initial DFT run.
		"""

		#generating wannier90.win
		TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
		TB.Compute_cor_idx(p['cor_at'],p['cor_orb'])
		print((TB.TB_orbs))
		if list(pV.keys()).count('NBANDS='):
			self.DFT.NBANDS = pV['NBANDS='][0]
		self.DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1])

	def vasp_run(self,dir):
		"""
		This method runs the inital VASP calculation.
		"""

		#initial VASP run
		print('\nRunning VASP in %s'%dir)
		cmd = 'cd '+dir+ ' && '+ self.para_com+' '+self.vasp_exec #+ " > dft.out 2> dft.error"
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if err:
			print('DFT calculation failed! Check dft.error for details.\n')
			errdir = dir+os.sep+'dft.error'
			f = open(errdir,'wb')
			f.write(err)
			f.close()
			sys.exit()
		else:
			print('DFT calculation complete.\n')	
			outdir = dir+os.sep+'dft.out'
			f = open(outdir,'wb')
			f.write(out)
			f.close()

	def vasp_convergence(self):

		"""
		This method checks for convergence inside the DFT_relax directory 
		and copies CONTCAR as POSCAR to root directory. Otherwise it runs vasp for convergence.
		If you want better convergence remember to copy an updated INCAR in the DFT_relax directory.
		"""

		#First check if  converged	
		if os.path.exists('./DFT_relax'):
			if os.path.exists('./DFT_relax/vasprun.xml'):
				vasprun = vasp.Vasprun('./DFT_relax/vasprun.xml')

				if vasprun.converged_ionic == True:
					print('Ionic convergence reached. Copying CONTCAR to root directory.')
					copyfile('./DFT_relax/CONTCAR','./POSCAR')

				else:	
					print('Convergence not reached. Recalculate with higher convergence parameters.')
					sys.exit()

			else:
				print('DFT_relax directory exists but vasprun.xml does not.')	
				sys.exit()
		else:
			print('DFT_relax directory does not exist.')
			sys.exit()

	def update_win(self):
		"""
		This updates the wannier90.win file with the number of bands and fermi energy 
		from the initial DFT calculation.
		"""

		#Updating wannier90.win with the number of DFT bands
		self.DFT.Read_NBANDS()
		self.DFT.Read_EFERMI()
		self.DFT.Update_win(self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1])


	def run_wan90(self):
		"""
		Running wannier90.x to generate .chk and .eig files.
		"""

		print('Running wannier90...')
		cmd = "wannier90.x wannier90"
		out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if err:
			print('wannier90 calculation failed!\n')
			print(err)
			sys.exit()
		else:
			print('wannier90 calculation complete.\n')	
			print(out)




	def gen_sig(self):
		"""
		This method generates the initial self energy file sig.inp.
		"""
		cmd = "sigzero.py"
		out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
		if err:
			print(err)
			sys.exit()
		else:
			print('Initial self-energy file generated.\n')


	def copy_files(self):	
		"""
		This creates a directory DMFT in the root directory 
		and copies all the necessary files.
		"""	

		#creating directory for DMFT
		if os.path.exists("DMFT"):
			shutil.rmtree("DMFT")
			os.makedirs("DMFT")
		else:	
			os.makedirs("DMFT")

		#copy INPUT.py to DMFT directory
		copyfile("INPUT.py","./DMFT/INPUT.py")

		#copying files into DMFT directory
		cmd = "cd ./DMFT && Copy_input_bands.py ../"
		out, err = subprocess.Popen(cmd, shell=True).communicate()
		if err:
			print('File copy failed!\n')
			print(err)
			sys.exit()
		else:
			print(out)
		print('\nDMFT initialization complete. Ready to run RUNDMFT.py.\n')

if __name__ == "__main__":

	#top level parser
	print('\n------------------------------- \n| DMFTwDFT Initialization tool |\n-------------------------------\n')
	des = 'This tool initiates the DMFT calculation by running DFT and wannier90.' 
	parser = argparse.ArgumentParser(description=des,formatter_class=RawTextHelpFormatter)
	

	#parser for dft  
	parser.add_argument('-dft',default='vasp', type=str, help='Choice of DFT code for the DMFT calculation. Currently available: \n [vasp,siesta]')
	parser.add_argument('-dftexec',default='vasp_std', type=str, help='The name of the chosed dft executable. \n E.g. vasp_std ')
	parser.add_argument('-relax',action='store_true', help='Flag to check for convergence. Program exits if not converged.')
	args = parser.parse_args()
	Initialize(args)
