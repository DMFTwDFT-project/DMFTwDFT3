#!/usr/bin/env python2
from scipy import *
import sys
import re, os, glob, shutil
from os.path import getsize


if __name__=='__main__':

   cpdr = sys.argv[1]
   post=None
   if (len(sys.argv)>2):
      post = sys.argv[2]

   if not os.path.exists(cpdr):
      print 'The directory '+cpdr+' does not exists!'
      sys.exit(1)

   if post=='-dos' or post=='dos':
      dosfiles=['wannier90.chk','wannier90.eig','DMFT_mu.out','dmft_params.dat','INPUT.py','para_com.dat']
      for files in dosfiles:
         if os.path.exists(cpdr+'/'+files):
            print "Copying dos file "+files+" to the current directory"
            shutil.copy2(cpdr+'/'+files, '.')
         else:
            print files+" must exist in a "+cpdr+" directory! Exiting!"; sys.exit(1)
      execfile('INPUT.py')
      if os.path.exists(p['path_bin']+'/Interpol_sig_real.py'):
         print "Copying dos file Interpol_sig_real.py to the current directory"
         shutil.copy2(p['path_bin']+'/Interpol_sig_real.py', '.')
      else:
         print "Interpol_sig_real.py must be copied from bin directory!"

   else:
      DFTfiles={'VASP':['OUTCAR','OSZICAR','POSCAR','POTCAR','KPOINTS','INCAR','WAVECAR','para_com.dat'],'QE':['XXX']}
      Ldft=False
      for dft in DFTfiles:
         outfile=DFTfiles[dft][0]
         if os.path.exists(cpdr+'/'+outfile): print dft+' results has been found in a '+cpdr+' directory!'; Ldft=True; break;
         if Ldft==False: print "No DFT results has been found in a "+cpdr+" directory!"#; break #sys.exit(1)

      if Ldft==True:
         for i,files in enumerate(DFTfiles[dft][:]):
            #copy files
            if os.path.exists(cpdr+'/'+files):
               print "Copying DFT file "+files+" to the current directory"
               shutil.copy2(cpdr+'/'+files, '.')
            else:
               if i<2:
                  print files+" must exist in a "+cpdr+" directory! Exiting!"; sys.exit(1)
               else:
                  print files+" does not exist in a "+cpdr+" directory!"
                  print files+" will be needed for charge updates!"
         if os.path.exists(cpdr+'/DFT_mu.out'):
            print "Copying DFT file DFT_mu.out to the current directory"
            shutil.copy2(cpdr+'/DFT_mu.out', '.')
         else:
            print "Making DFT_mu.out file"
            fi=open(cpdr+'/'+outfile,'r')
            EFermi=0
            for line in fi:
               if re.search('Fermi energy',line) or re.search('E-fermi',line):
                  line_fermi=line
            #print line_fermi
            val=re.search(r'(\-?\d+\.?\d*)',line_fermi)
            #print val
            EFermi=float(val.group(1))
            savetxt('./DFT_mu.out',array([EFermi]))


      Wannierfiles=['wannier90.chk','wannier90.eig','wannier90.win','wannier90.amn']
      for i,files in enumerate(Wannierfiles):
         #copy files
         if os.path.exists(cpdr+'/'+files):
            print "Copying Wannier file "+files+" to the current directory"
            shutil.copy2(cpdr+'/'+files, '.')
         else:
            print files+" must exist in a "+cpdr+" directory! Exiting!"; sys.exit(1)
            #if i<2:
            #   print files+" must exist in a "+cpdr+" directory! Exiting!"; sys.exit(1)
            #else:
            #   print files+" does not exist in a "+cpdr+" directory!"
            #   print files+" will be needed for charge update!"
               

      DMFTfiles=['sig.inp','DMFT_mu.out']
      if os.path.exists(cpdr+'/sig.inp'):
         print "Copying DMFT file sig.inp to the current directory"
         shutil.copy2(cpdr+'/sig.inp', '.')
      else:
         print "sig.inp file does not exist! Must be generated using sigzero.py"
      if os.path.exists(cpdr+'/DMFT_mu.out'):
         print "Copying DMFT file DMFT_mu.out to the current directory"
         shutil.copy2(cpdr+'/DMFT_mu.out', '.')
      else:
         print "DMFT_mu.out file does not exist! Copying from DFT_mu.out"
         shutil.copy2('./DFT_mu.out', './DMFT_mu.out')


