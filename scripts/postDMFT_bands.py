#!/usr/bin/env python
import sys, subprocess, os
import numpy as np
from scipy import *
import copy, Fileio, re
from scipy import interpolate
import shutil
import glob
from INPUT import *
import argparse

class PostProcess:
	"""DMFTwDFT PostProcess

	This class contains methods to perform post processing of the DMFT calculations. 

	inputs:
	emin
	emax
	rom : omega points
	broaden
	siglistindx: How many self-energies?

	"""

	def __init__(self,emin,emax,rom,broaden,siglistindx):
		"""
		Initializes the following:
		"""

		#mpirun
		if os.path.exists("para_com.dat"):
    		fipa=open('para_com.dat','r')
    		self.para_com=str(fipa.readline())[:-1]
    		fipa.close()
		else:
    		self.para_com=""


   		self.emin = None
   		self.emax = None
   		self.rom = None
   		self.broaden = None 	

		print('\n')
		print('-------------------------------')
		print('| DMFTwDFT post-processing scheme |')
		print('-------------------------------\n')

		#dmft_bin
		self.path_bin=p['path_bin']


	def anal_cont(self):
		"""
		This method performs the analytic continuation.
		"""

		#creating directory for ac
		if os.path.exists("ac"):
			shutil.rmtree("ac")
			os.makedirs("ac")
		else:	
			os.makedirs("ac")			

		#copying the last few self-energies from the DMFT run in the directory above
		siglist = sorted(glob.glob("sig.inp.*"),key=os.path.getmtime)[-self.siglistindx:]
		for file in siglist:
			shutil.copy(file,'ac')

		#averaging sef energies
		print('Averaging self-energies from: ')
		print('%s'siglist)
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
		    shutil.copyfile('maxent_params.dat','./ac/maxent_params.dat')
		else:
		    src=path_bin+ os.sep+"maxent_params.dat"
		    shutil.copyfile(src,"./ac/maxent_params.dat")

		#Analytic continuation
		print('Running analytic continuation...')
		cmd = "cd ac && maxent_run.py sig.inpx"+" > ac.out 2> ac.error"
		out, err = subprocess.Popen(cmd, shell=True,stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if os.path.exists("./ac/Sig.out"):
		    print('Analytic continuation complete.\n')  
		else:
		    print('Analytic continuation failed! Check ac.error for details.\n') 
		    sys.exit()
  


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='DMFTwDFT post-processing tool')
	parser.add_argument('-emin',default=-7.0, type=float, help='Minimum value for interpolation')
	parser.add_argument('emax',default=3.0, type=float, help='Maximum value for interpolation')
	parser.add_argument('rom',default=3000, type=float, help='Matsubara Frequency (omega) points')
	parser.add_argument('broaden',default=0.03, type=float, help='Broadening')
	parser.add_argument('siglistindx',default=2, type=int, help='How many last self energy files to average?')
	args = parser.parse_args()

	PostProcess(args.emin,args.emax,args.rom,args.broaden,args.siglistindx)
