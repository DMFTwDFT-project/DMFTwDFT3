#!/usr/bin/env python

from scipy import *
from scipy import interpolate
import sys
sys.path.insert(0,"/home/uthpala/Documents/Research/CondensedMatterResearch/projects/DFTDMFT/bin/")
import Struct, copy, Fileio


if __name__=='__main__':

   emin=-5.0;emax=5.0;rom=100
   execfile('INPUT.py')
   TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
   cor_at=p['cor_at']; cor_orb=p['cor_orb']
   TB.Compute_cor_idx(cor_at,cor_orb)
   ommesh=linspace(emin,emax,rom)
   Sig_tot=zeros((TB.ncor_orb,rom),dtype=complex)
   for i,ats in enumerate(cor_at):
      (om,Sig,TrSigmaG,Epot,nf_q,mom)=Fileio.Read_complex_Data('Sig'+str(i+1)+'.out')
      newSig=zeros((len(Sig),rom),dtype=complex)
      for ii in range(len(Sig)):
         SigSpline = interpolate.splrep(om, Sig[ii].real, k=1, s=0)
         newSig[ii,:] += interpolate.splev(ommesh, SigSpline)
         SigSpline = interpolate.splrep(om, Sig[ii].imag, k=1, s=0)
         newSig[ii,:] += 1j*interpolate.splev(ommesh, SigSpline)
      for at in ats:
         for ii,orbs in enumerate(cor_orb[i]):
            for orb in orbs: 
               idx=TB.idx[at][orb]
               Sig_tot[idx,:]=copy.deepcopy(newSig[ii,:])
   Fileio.Print_complex_multilines(Sig_tot,ommesh,'SigMoo_real.out')
