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
	

	if the -relax flag is used then the program will check if convergence is reached and rerun
	if necessary. Remember to put an updated version of the DFT inputs for higher convergence 
	inside DFT_relax.

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
		
		if args.dmft:
			self.run_dmft(args.force)

		if args.hf:
			self.run_hf(args.force)

	def gen_win(self):
		"""
		This method generates wannier90.win for initial DFT run.
		"""
		#wannier mesh tolerance 
		kmesh_tol = 0.000001

		#generating wannier90.win
		TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
		TB.Compute_cor_idx(p['cor_at'],p['cor_orb'])
		print((TB.TB_orbs))
		if list(pV.keys()).count('NBANDS='):
			self.DFT.NBANDS = pV['NBANDS='][0]
		self.DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],self.DFT.NBANDS,self.DFT.EFERMI+p['ewin'][0],self.DFT.EFERMI+p['ewin'][1],kmesh_tol)

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
				try:
					vasprun = vasp.Vasprun('./DFT_relax/vasprun.xml')

					if vasprun.converged_ionic == True:
						print('Ionic convergence reached. Copying CONTCAR to root directory.')
						copyfile('./DFT_relax/CONTCAR','./POSCAR')

					else:	
						print('Convergence not reached. Recalculating with higher convergence parameters.')
						self.vasp_run('./DFT_relax')	

						try:
							vasprun = vasp.Vasprun('./DFT_relax/vasprun.xml')
							if vasprun.converged_ionic == True:
								print('Ionic convergence reached. Copying CONTCAR to root directory.')
								copyfile('./DFT_relax/CONTCAR','./POSCAR')
							else:
								print('Convergence not reached. Update convergence parameters.')
								sys.exit()
						except:
							print('vasprun.xml is incomplete. This is a result of an incomplete VASP calculation.')	
							sys.exit()
				except:
					print('vasprun.xml is incomplete. This is a result of an incomplete VASP calculation.')		
					print('Restarting VASP calculation...')	
					self.vasp_run('./DFT_relax')
					try:
						vasprun = vasp.Vasprun('./DFT_relax/vasprun.xml')	
						if vasprun.converged_ionic == True:
							print('Ionic convergence reached. Copying CONTCAR to root directory.')
							copyfile('./DFT_relax/CONTCAR','./POSCAR')
						else:
							print('Convergence not reached. Update convergence parameters.')
							sys.exit()
					except:
						print('vasprun.xml is incomplete. This is a result of an incomplete VASP calculation.')		
						sys.exit()
			else:
				print('DFT_relax directory exists but vasprun.xml does not. Running VASP now.')	
				self.vasp_run('./DFT_relax')


				try:
					vasprun = vasp.Vasprun('./DFT_relax/vasprun.xml')	
					if vasprun.converged_ionic == True:
						print('Ionic convergence reached. Copying CONTCAR to root directory.')
						copyfile('./DFT_relax/CONTCAR','./POSCAR')
					else:
						print('Convergence not reached. Update convergence parameters.')
						sys.exit()
				except:
					print('vasprun.xml is incomplete. This is a result of an incomplete VASP calculation.')		
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
			if os.path.exists("DMFT/imp.0/"):
				shutil.rmtree("DMFT/imp.0/")
				#os.makedirs("DMFT")
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
		print('\nDMFT initialization complete. Ready to run DMFT calculation.\n')

	def copy_files_hf(self):	
		"""
		This creates a directory HF in the root directory 
		and copies all the necessary files.
		"""	

		#creating directory for HF 
		if os.path.exists("HF"):
			if os.path.exists("HF/imp.0/"):
				shutil.rmtree("HF/imp.0/")
				#os.makedirs("DMFT")
		else:	
			os.makedirs("HF")

		#copy INPUT.py to DMFT directory
		copyfile("INPUT.py","./HF/INPUT.py")

		#copying files into DMFT directory
		cmd = "cd ./HF && Copy_input_bands.py ../"
		out, err = subprocess.Popen(cmd, shell=True).communicate()
		if err:
			print('File copy failed!\n')
			print(err)
			sys.exit()
		else:
			print(out)
		print('\nHF initialization complete. Ready to run HF calculation.\n')


	def run_dmft(self,force):
		"""
		This first checks if there is a previous DMFT calculation and runs
		DMFT only if that run is incomplete unless forced.
		"""

		#Checking for previous DMFT run in the directory
		pathstr ='./DMFT'+os.sep+'INFO_TIME'

		if os.path.exists(pathstr):
			fi=open(pathstr,'r')
			done_word=fi.readlines()[-1]
			fi.close()

			if done_word.split()[0] == 'Calculation':
				print('Existing DMFT calculation is complete.')
				if force:
					#forcing DMFT calculation	
					print('-force flag enabled. Restarting DMFT...')
					self.vasp_run(self.dir)
					self.update_win()
					self.run_wan90()
					self.copy_files()
					print('Running DMFT...')
					cmd = 'cd '+'DMFT'+ ' && '+ 'RUNDMFT.py' #+ " > dft.out 2> dft.error"
					out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					if err:
						print('DMFT calculation failed! Check dmft.error for details.\n')
						errdir = './DMFT'+os.sep+'dmft.error'
						f = open(errdir,'wb')
						f.write(err)
						f.close()
						sys.exit()
					else:
						print('DMFT calculation complete.\n')	
						outdir ='./DMFT'+os.sep+'dmft.out'
						f = open(outdir,'wb')
						f.write(out)
						f.close()

				else:
					#exit when exiting DMFT calculation is complete.
					print('-force flag disabled. Exiting. ')
					sys.exit()


			else:
				#Incomplete DMFT calculation.
				print('DMFT calculation incomplete.')
				self.vasp_run(self.dir)
				self.update_win()
				self.run_wan90()
				self.copy_files()
				print('Running DMFT...')
				cmd = 'cd '+'DMFT'+ ' && '+ 'RUNDMFT.py' #+ " > dft.out 2> dft.error"
				out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				if err:
					print('DMFT calculation failed! Check dmft.error for details.\n')
					errdir = './DMFT'+os.sep+'dmft.error'
					f = open(errdir,'wb')
					f.write(err)
					f.close()
					sys.exit()
				else:
					print('DMFT calculation complete.\n')	
					outdir ='./DMFT'+os.sep+'dmft.out'
					f = open(outdir,'wb')
					f.write(out)
					f.close()

		else:
			#no DMFT/INFO_TIME found
			self.vasp_run(self.dir)
			self.update_win()
			self.run_wan90()
			self.copy_files()
			print('Running DMFT...')
			cmd = 'cd '+'DMFT'+ ' && '+ 'RUNDMFT.py' #+ " > dft.out 2> dft.error"
			out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			if err:
				print('DMFT calculation failed! Check dmft.error for details.\n')
				errdir = './DMFT'+os.sep+'dmft.error'
				f = open(errdir,'wb')
				f.write(err)
				f.close()
				sys.exit()
			else:
				print('DMFT calculation complete.\n')	
				outdir ='./DMFT'+os.sep+'dmft.out'
				f = open(outdir,'wb')
				f.write(out)
				f.close()

	def run_hf(self,force):
		"""
		This first checks if there is a previous HF calculation and runs
		HF only if that run is incomplete unless forced.
		"""

		#Checking for previous HF run in the directory
		pathstr ='./HF'+os.sep+'INFO_TIME'

		if os.path.exists(pathstr):
			fi=open(pathstr,'r')
			done_word=fi.readlines()[-1]
			fi.close()

			if done_word.split()[0] == 'Calculation':
				print('Existing HF calculation is complete.')
				if force:
					#forcing HF calculation	
					print('-force flag enabled. Restarting HF...')
					self.vasp_run(self.dir)
					self.update_win()
					self.run_wan90()
					self.copy_files_hf()
					print('Running HF...')
					cmd = 'cd '+'HF'+ ' && '+ 'RUNDMFT_HF.py' #+ " > dft.out 2> dft.error"
					out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					if err:
						print('HF calculation failed! Check hf.error for details.\n')
						errdir = './HF'+os.sep+'hf.error'
						f = open(errdir,'wb')
						f.write(err)
						f.close()
						sys.exit()
					else:
						print('HF calculation complete.\n')	
						outdir ='./HF'+os.sep+'hf.out'
						f = open(outdir,'wb')
						f.write(out)
						f.close()

				else:
					#exit when exiting HF calculation is complete.
					print('-force flag disabled. Exiting. ')
					sys.exit()


			else:
				#Incomplete HF calculation.
				print('HF calculation incomplete.')
				self.vasp_run(self.dir)
				self.update_win()
				self.run_wan90()
				self.copy_files_hf()
				print('Running HF...')
				cmd = 'cd '+'HF'+ ' && '+ 'RUNDMFT_HF.py' #+ " > dft.out 2> dft.error"
				out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				if err:
					print('HF calculation failed! Check hf.error for details.\n')
					errdir = './HF'+os.sep+'hf.error'
					f = open(errdir,'wb')
					f.write(err)
					f.close()
					sys.exit()
				else:
					print('HF calculation complete.\n')	
					outdir ='./HF'+os.sep+'hf.out'
					f = open(outdir,'wb')
					f.write(out)
					f.close()

		else:
			#no HF/INFO_TIME found
			self.vasp_run(self.dir)
			self.update_win()
			self.run_wan90()
			self.copy_files_hf()
			print('Running HF...')
			cmd = 'cd '+'HF'+ ' && '+ 'RUNDMFT_HF.py' #+ " > dft.out 2> dft.error"
			out, err = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			if err:
				print('HF calculation failed! Check hf.error for details.\n')
				errdir = './HF'+os.sep+'hf.error'
				f = open(errdir,'wb')
				f.write(err)
				f.close()
				sys.exit()
			else:
				print('HF calculation complete.\n')	
				outdir ='./HF'+os.sep+'hf.out'
				f = open(outdir,'wb')
				f.write(out)
				f.close()			


if __name__ == "__main__":

	#top level parser
	print('\n------------------------------- \n| DMFTwDFT Initialization tool |\n-------------------------------\n')
	des = 'This tool initializes the DMFT calculation by running DFT and wannier90 and the runs DMFT.' 
	parser = argparse.ArgumentParser(description=des,formatter_class=RawTextHelpFormatter)
	

	#parser for dft  
	parser.add_argument('-dft',default='vasp', type=str, help='Choice of DFT code for the DMFT calculation. Currently available: \n [vasp,siesta]')
	parser.add_argument('-dftexec',default='vasp_std', type=str, help='The name of the chosed dft executable. \n E.g. vasp_std ')
	parser.add_argument('-relax',action='store_true', help='Flag to check for DFT convergence. Program exits if not converged.')
	parser.add_argument('-dmft',action='store_true',help='Flag to run DMFT. Checks for a previous DMFT calculation and runs only if it is incomplete.')
	parser.add_argument('-hf',action='store_true',help='Flag to perform Hartree-Fock calculation to the correlated orbitals.')
	parser.add_argument('-force',action='store_true',help='Flag to force DMFT or HF calculation even if a previous calculation has been completed.')
	args = parser.parse_args()
	Initialize(args)
