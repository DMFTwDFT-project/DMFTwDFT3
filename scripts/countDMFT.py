#!/usr/bin/env python

import sys,os
import argparse

def count_complete(args):
	"""
	This methods checks to see if the DMFT calculation is done
	by checking if "Done" is printed in INFO_TIME. It also checks
	for completeness in post-processing calculations. 

	Parameters
	----------
	path : str
		The Path to the DMFT/HF directory. 
		Could be the root directory of many calculations as well.

	type : str
		The type of calculation. Either DMFT or HF.

	post : str
		Checks for completeness of post-processing calculations. 
		ac, dos or bands. 


	"""

	filedic={
		'ac' : 'Sig.out',
		'dos' : 'DMFT-PDOS.png',
		'plainbands' : 'A_k.eps',
		'partialbands' : 'A_k_partial.eps',
	}

	done_counter = 0

	try:
		pathlist = sorted([int(d) for d in os.listdir(args.path) if os.path.isdir(d)])
		print(pathlist)
	except:
		pathlist = [os.getcwd()]

	for path in pathlist:

		pathstr = str(path)+os.sep+args.type.upper()+os.sep+'INFO_TIME'

		if os.path.exists(pathstr):
			fi=open(pathstr,'r')
			done_word=fi.readlines()[-1]
			fi.close()

			if done_word.split()[0] == 'Calculation':
				done_counter += 1
				print('%s calculation complete at %s' %(args.type.upper(),path))

				if args.post:
					for i in args.post:
						if i == 'plainbands' or i == 'partialbands':
							ii = 'bands'
						else:
							ii = i
						postpathstr = str(path)+os.sep+args.type.upper()+os.sep+ii+os.sep+filedic[i]
						if os.path.exists(postpathstr):
							print('%s complete.' %i)	
						else:
							print('%s incomplete.' %i)	
					print('\n')	


			else:
				print('%s calculation incomplete at %s' %(args.type.upper(),path))


		else:
			print('INFO_TIME does not exist at %s' %path)

	
	print('%d %s calculations have been completed.'%(done_counter,args.type.upper()))			

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='This script checks to see if the DMFT/HF calculation is complete.')
	parser.add_argument('path', type=str,default='./', help='Path to DMFT/HF directory')
	parser.add_argument('-type',type=str,default='dmft',help='DMFT or HF',choices=['dmft','hf'])
	parser.add_argument('-post',type=str,default=None,help='Check for post-processing completeness.',choices=['ac','dos','plainbands','partialbands'],nargs='+')
	args = parser.parse_args()
	count_complete(args)

		
