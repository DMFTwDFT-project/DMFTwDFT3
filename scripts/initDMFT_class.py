#!/usr/bin/env python2
import sys, subprocess, os
import numpy as np
import shutil
from shutil import copyfile
import VASP
import Struct
from INPUT import *
import argparse
from argparse import RawTextHelpFormatter


class Initialize():
	"""DMFTwDFT Initialization

	This class contains methods to run the initial DFT and wannier90 calculation
	to generate inputs for the DMFT runs.
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
			self.DFT=VASP.VASP_class()

		
		self.gen_win()	
		self.gen_sig()
		self.vasp_run()
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
		print(TB.TB_orbs)
		if pV.keys().count('NBANDS='):
			self.DFT.NBANDS = pV['NBANDS='][0]
		self.DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1])

	def vasp_run(self):
		"""
		This method runs the inital VASP calculation.
		"""

		#initial VASP run
		print('\nRunning initial VASP...')
		cmd = self.para_com+" "+self.vasp_exec #+ " > dft.out 2> dft.error"
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
		cmd = "cd ./DMFT && Copy_input.py ../"
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
	
	args = parser.parse_args()
	Initialize(args)
