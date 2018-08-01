#!/usr/bin/env python


import sys
sys.path.insert(0, "/home/uthpala/Documents/Research/CondensedMatterResearch/projects/DFTDMFT/bin/")

#from data import *
import sys, subprocess
import re, os, glob, shutil, socket, time
from os.path import getsize
import Struct,Fileio,VASP, DMFT_MOD, IMP_SOLVER
from scipy import *
import Read_Prob, scipy.interpolate

#edits
#commented out if loop from 311-341

#parallel run
para_com=""

#vasp executable
vasp_exec="/home/uthpala/Documents/Research/CondensedMatterResearch/projects/DFTDMFT/bin/vasp_dmft"

###########################################################################
#### This program executes VASP+DMFT using CTQMC impurity solver ##########
#### This part controls the overall self-consistent steps        ##########
############################################  Dr. Hyowon Park    ##########

def now(): return time.strftime(' at %H:%M:%S')

def Print_dm_out(DMFT,TB,cor_at):

      ####### Print final results ##############
   fi=open('dm.out','w')
   print >>fi, 'mu=',DMFT.mu
   print >>fi, 'E_dc=',
   for i in range(len(cor_at)):
      print >>fi, '%5f ' %(DMFT.VDC[i]),
   print >>fi, ''
   Nd_eg=zeros(len(cor_at),dtype=float); Nd_t2g=zeros(len(cor_at),dtype=float)
   print >>fi, "Nd_eg=",
   for i,ats in enumerate(cor_at):
      for at in ats:
         Nd_eg[i]+=float(DMFT.dm[TB.idx[at]['d_z2']])+float(DMFT.dm[TB.idx[at]['d_x2y2']])
      print >>fi, Nd_eg[i]/len(ats),
   print >>fi, ''
   print >>fi, "Nd_t2g=",
   for i,ats in enumerate(cor_at):
      for at in ats:
         Nd_t2g[i]+=float(DMFT.dm[TB.idx[at]['d_xz']])+float(DMFT.dm[TB.idx[at]['d_yz']])+float(DMFT.dm[TB.idx[at]['d_xy']])
      print >>fi, Nd_t2g[i]/len(ats),
   print >>fi, ''
   print >>fi, "Nd=",
   for i,ats in enumerate(cor_at):
      print >>fi, (Nd_eg[i]+Nd_t2g[i])/len(ats),
   print >>fi, ''
   for i in range(len(DMFT.dm)): print >>fi, '%8.6f ' %(DMFT.dm[i]),
   print >>fi, ''
   for i in range(len(DMFT.MOM)): print >>fi, '%8.6f ' %(DMFT.MOM[i]),
   print >>fi, ''
   fi.close()
#   fi=open('sigoo.out','w')
#   print >>fi, 'mu=',mu
#   print >>fi, 'E_dc=',
#   for i in range(len(cor_at)):
#      print >>fi, '%5f ' %(E_dc[i]),
#   print >>fi, ''
#   for i in range(len(sigi)): print >>fi, '%5f ' %(sigi[i,-1].real),
#   print >>fi, ''
#   fi.close()
#

def CreateINCAR(params_vasp):
   " Creates input file (INCAR) for vasp"
   f = open('INCAR', 'w')
   for p in params_vasp:
       print >> f, p, '  ', params_vasp[p][0], '\t', params_vasp[p][1]
   f.close()   


#class Parameters:
#"""Class to store various paramers to control the flow of the program
#   Input:
#      f_sparams    -- filename of the parameters read only at the beginning
#      sparams      -- a dictionary of pairs {name of the variable: default value} which are looked for in the f_sparams file
#      f_params     -- filename of the parameters which can be refreshed several times by refresh member function
#      params       -- a dictionary of {variable: default value} pairs which can be updated at runtime
#   Output:
#      None
#   Variables which become members of the class are those listed in params and sparams dictionaries
#"""
#   def __init__(self, ):


if __name__ == '__main__':
   

   execfile('INPUT.py') # Read input file
   main_out=open('TIMEINFO','w')
   main_iter=open('ENERGY','w')
   DMFT_iter=open('ITERINFO','w')
   DM_iter=open('DMINFO','w')
   E_iter=open('ENERGY_AVG','w')


   if os.path.exists("para_com.dat"):
      fipa=open('para_com.dat','r')
      para_com=str(fipa.readline())[:-1]
      fipa.close()
   else:
      para_com=para_com

 ############ Initial Preparation ########################

   if p['Nit']>0 and p['Niter']>1: main_out.write( '-----------Charge Self-consistent DFT+DMFT calculations-----------' )
   if p['Nit']>0 and p['Niter']==1: main_out.write( '----------Non-charge Self-consistent DFT+DMFT calculations-----------' )
   main_out.write('\n')
   main_out.flush()

   main_out.write( 'Caculation Starts'+now() )
   main_out.write('\n')
   main_out.flush()

   DMFT_iter.write( '%6s %8s %8s %8s %8s %8s %8s %8s' % ('mu','TOTN','Nd[0]','Nd[-1]','EKIN','EPOT','Sig[0]','Sig[-1]') )
   #for iat in p['cor_at']:   # Over all inequivalent impurity problems
   #    DMFT_iter.write( '%12s %12s %12s %12s %12s ' % ('n_latt', 'n_imp','Sigoo[0]','Sigoo[-1]','VDC') )
   DMFT_iter.write('\n')
   DMFT_iter.flush()

   DM_iter.write( '%10s' % ('Occupancy') )
   #for iat in p['cor_at']:   # Over all inequivalent impurity problems
   #    DMFT_iter.write( '%12s %12s %12s %12s %12s ' % ('n_latt', 'n_imp','Sigoo[0]','Sigoo[-1]','VDC') )
   DM_iter.write('\n')
   DM_iter.flush()

   main_iter.write( '%3s %5s %12s %12s %12s %12s' % ('N','NDMFT','TOTENERGY','EDIFF','TOTENERGY2','EDIFF') )
   main_iter.write('\n')
   main_iter.flush()

   #### Read POSCAR ########
   TB=Struct.TBstructure('POSCAR',p['atomnames'],p['orbs'])
   cor_at=p['cor_at']; cor_orb=p['cor_orb']
   TB.Compute_cor_idx(cor_at,cor_orb)
   print TB.TB_orbs
   if TB.LHF=='.TRUE.': p['Nd_qmc']=0
   U=p['U'];J=p['J']
   T_high=p['T_high']
   if T_high<1.0/pC['beta'][0]: T_high=1.0/pC['beta'][0]

   print p['noms']/(T_high*pC['beta'][0])
   
   noms_high=int(p['noms']/(T_high*pC['beta'][0]))

   DMFT=DMFT_MOD.DMFT_class(p,pC,TB)

   DFT=VASP.VASP_class()
   ETOT_old=0.0;ETOT2_old=0.0;ETOT=0.0;ETOT2=0.0
   CHGDIFF=0.;CHGDIFF2=0.

   for itt in range(p['Niter']):
      main_out.write( '--- Starting Charge loop '+str(itt+1)+now()+'---' )
      main_out.write('\n')
      main_out.flush()

      if itt<p['Niter']-p['Nforce']: Fileio.Creat_INPUT(p,pC,TB,T_high,noms_high,'.FALSE.')
      else: Fileio.Creat_INPUT(p,pC,TB,T_high,noms_high,'.TRUE.')

      if p['Nrelax']>p['Nforce']: print "Nrelax should not be larger than Nforce"; exit()
      if itt>=p['Niter']-p['Nrelax']:
         TB.Update_POSCAR('POSCAR')
         pV['NSW= ']=[3,'# NSW']
         pV['POTIM= ']=[0.1,'# NSW']
         pV['IBRION= ']=[1,'# IBRION']
         pV['ISIF= ']=[2,'# ISIF']
         pV['EDIFFG= ']=[-0.01,'# EDIFFG']
      CreateINCAR(pV)
      if os.path.exists('CHGCAR') and itt==0:
         f = open('INCAR', 'a')
         print >> f, "ICHARG= 11"
         f.close()


      main_out.write( '--- Running vaspCHG '+now()+'---' )
      main_out.write('\n')
      main_out.flush()

      #added on July 18,2018 because LDMFT doesn't work for first DFT run
      if itt>0:
         f = open('INCAR', 'a')
         print >> f, "LDMFT= .TRUE."
         print >> f, "NELM= 20"
         f.close()

      if itt==0:
         if pV.keys().count('NBANDS='): DFT.NBANDS=pV['NBANDS='][0]
      DFT.Create_win(TB,p['atomnames'],p['orbs'],p['L_rot'],DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1],0.000001)
      cmd = vasp_exec+" > vasp.out 2> vasp.error"
      #cmd = p['path_bin']+"vaspCHG > vasp.out 2> vasp.error"
      out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
      print os.popen("cp vasp.out vasp.out."+str(itt)).read()
      print os.popen("cp OUTCAR OUTCAR."+str(itt)).read()
      print os.popen("cp OSZICAR OSZICAR."+str(itt)).read()
      print os.popen("cp CHGCAR CHGCAR."+str(itt)).read()
      if itt>=p['Niter']-p['Nrelax']: 
         print os.popen("cp POSCAR POSCAR."+str(itt)).read()
         print os.popen("cp CONTCAR POSCAR").read()

      if itt>0:
         DFT.Read_NELECT()
         CHGDIFF=DFT.Diff_CHGCAR('CHGCAR.'+str(itt-1),'CHGCAR.'+str(itt))
      DFT.Read_NBANDS()
      DFT.Read_EFERMI()
      DFT.Update_win(DFT.NBANDS,DFT.EFERMI+p['ewin'][0],DFT.EFERMI+p['ewin'][1])
      if itt==0:
         fi=open('DMFT_mu.out','w')
         print >>fi, DFT.EFERMI
         fi.close()

      if itt>0:# or (not os.path.exists('CHGCAR')):
         print os.popen("rm wannier90.chk").read()
         print os.popen("rm wannier90.chk.fmt").read()
#         print os.popen("rm wannier90.eig").read()

      #if itt<p['Niter']-1:
      main_out.write( '-------------- Running wannier 90 '+str(itt+1)+'----------------' )
      main_out.write('\n')
      main_out.flush()
      cmd = p['path_bin']+"wannier90.x wannier90"
      out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
      print out #, err


      for it in range(p['Nit']):
         main_out.write( '--- Starting DMFT loop '+str(it+1)+now()+'---' )
         main_out.write('\n')
         main_out.flush()
         
         if itt==0 and it==0:
            for i,ats in enumerate(cor_at):
               if p['nspin']==2: pD['para=']=0 
               else: pD['para=']=1
               if p['orbs'][i]=='f': pD['l=']=3
               elif p['orbs'][i]=='d': pD['l=']=2
               pD['J=']=float(J[i])
               pD['Eimp=']=zeros(p['nspin']*len(cor_orb[i]))#array(ed[i])-ed[i][0]+array(DMFT.sig_st[i])-DMFT.sig_st[i][0]
               IMP_SOLVER.Create_atomd(pD)
               IMP_SOLVER.Create_Trans(TB.ncor_orb_list[i],p['nspin'],ats[0],cor_orb[i],TB)
               cmd = 'cp Trans.dat Trans'+str(i+1)+'.dat'
               print os.popen(cmd).read()
               cmd = p['path_bin']+"atom_d.py 'atom_d.inp' 2> atom_d.error || { echo 'ATOMD run failed!'; exit 1; }"
               out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
               print out #, err
               cmd = 'cp UC.dat UC'+str(i+1)+'.dat'
               print os.popen(cmd).read()
            DMFT.Read_Sig(TB,p['nspin'])
            DMFT.Compute_HF(1,p,TB)
            DMFT.Cmp_Sig_highT(T_high,noms_high,p['nspin'])
            #DMFT.Mix_Sig(p['mix_sig'])
            if TB.LHF=='.FALSE.':
               Fileio.Print_complex_multilines(DMFT.Sig,DMFT.ommesh,'SigMoo.out')
               if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn,DMFT.ommesh,'SigMoo_dn.out')
               Fileio.Print_complex_multilines(DMFT.Sig_highT,DMFT.ommesh_highT,'SigMoo_highT.out')
               if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn_highT,DMFT.ommesh_highT,'SigMoo_dn_highT.out')

            if (os.path.exists('SigMdc.out')):
               pass
            else:
               fi=open('SigMdc.out','w')
               for data in DMFT.SigMdc:
                  print >>fi, data,
               print >>fi, ''
               fi.close()
            if p['nspin']==2: 
               if (os.path.exists('SigMdc_dn.out')):
                  pass
               else:
                  fi=open('SigMdc_dn.out','w')
                  for data in DMFT.SigMdc_dn:
                     print >>fi, data,
                  print >>fi, ''
                  fi.close()
   
   
         if TB.LHF=='.FALSE.':
            cmd = para_com+" "+p['path_bin']+"dmft_ksum_sp > ksum_output 2> ksum_error"
         else: 
            cmd = para_com+" "+p['path_bin']+"XHF.py > ksum_output 2> ksum_error"
         out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
         print out #, err
   
         if TB.LHF=='.FALSE.': DMFT_mu,ed=DMFT.Compute_Delta(DMFT.T,p['nspin'],cor_at,cor_orb,TB,p['nom'])
         
         if itt==0 and it==0:
            cmd = 'cp DMFT_mu.out VASP_mu.out'
            print os.popen(cmd).read()
            #fi=open('VASP_mu.out','w')
            #print >>fi, DMFT_mu 
            #fi.close()
   
         main_out.write( '--- DMFT Ksum is finished '+now()+'---' )
         main_out.write('\n')
         main_out.flush()
   
         DMFT.Update_Nlatt(TB,p['nspin'])
         if not (os.path.exists('Sig1.out')): 
            DMFT.Nd_imp=DMFT.Nd_latt
            DMFT.Compute_HF(p['Nd_qmc'],p,TB)
   
         #DMFT.Impurity_Solver
         if TB.LHF=='.FALSE.': IMP_SOLVER.RUN_CTQMC(p,pC,pD,it,itt,para_com,DMFT_mu,ed,DMFT.sig_st)
         main_out.write( '--- Impurity solver is finished '+now()+'---' )
         main_out.write('\n')
         main_out.flush()
    
         DMFT.Read_Sig(TB,p['nspin'])
         DMFT.Compute_HF(p['Nd_qmc'],p,TB)
         DMFT.Mix_Sig(p['nspin'],p['mix_sig'])
         DMFT.Cmp_Sig_highT(T_high,noms_high,p['nspin'])
         if TB.LHF=='.FALSE.':
            Fileio.Print_complex_multilines(DMFT.Sig,DMFT.ommesh,'SigMoo.out')
            if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn,DMFT.ommesh,'SigMoo_dn.out')
            Fileio.Print_complex_multilines(DMFT.Sig_highT,DMFT.ommesh_highT,'SigMoo_highT.out')
            if p['nspin']==2: Fileio.Print_complex_multilines(DMFT.Sig_dn_highT,DMFT.ommesh_highT,'SigMoo_dn_highT.out')
         fi=open('SigMdc.out','w')
         for data in DMFT.SigMdc:
            print >>fi, data,
         print >>fi, ''
         fi.close()
         if p['nspin']==2: 
            fi=open('SigMdc_dn.out','w')
            for data in DMFT.SigMdc_dn:
               print >>fi, data,
            print >>fi, ''
            fi.close()

   
         fiDMFT=open('ITERINFO','r')
         Eline=fiDMFT.readlines()[-1:]
         E_KIN=float(Eline[0].split()[4])
         E_POT=float(Eline[0].split()[5])
         E_POT2=0.0
         for i,ats in enumerate(cor_at):
            E_POT2+=len(ats)*DMFT.Eimp[i]
         fiDMFT.close()
         ECOR=E_KIN
         #if itt>0:
         #   fiEKIN=open('VASP.output','r')
         #   ECOR=float(fiEKIN.readline().split()[3])
         #   fiEKIN.close()

         DFT.Read_OSZICAR('OSZICAR.'+str(itt)) 
         #read VASP.log instead of VASP.output 
         if itt==0:
            ETOT=DFT.E+E_POT
            ETOT2=DFT.E+E_POT
         else:
            fiEKIN=open('VASP.log','r')
            lines=fiEKIN.readlines()
            EB=float(lines[-1].split()[0])
            ECOR=float(lines[-1].split()[1])
            fiEKIN.close()
            ETOT=DFT.E-EB+E_KIN+E_POT#+ECOR #-ECOR+E_KIN+E_POT
            ETOT2=DFT.E+E_POT         


         ETOT_old=ETOT; ETOT2_old=ETOT2
         VdcNd=0
         for i,ats in enumerate(cor_at):
            VdcNd+=len(ats)*DMFT.VDC[i]*DMFT.Nd_imp[i]
         #ETOT=DFT.E-ECOR+E_KIN+E_POT
         ETOT=DFT.E+E_POT
         #ETOT=DFT.E-ECOR+E_KIN+E_POT2+DMFT.EHF_latt-DMFT.EHF_imp-DMFT.EDC
         #ETOT2=DFT.E-ECOR+E_KIN+E_POT2-DMFT.EHF_imp+DMFT.EHF_latt-DMFT.EDC
         ETOT2=DFT.E-ECOR+E_KIN+E_POT2-DMFT.EDC_imp
         #ETOT2=DFT.E-ECOR+E_KIN+2*E_POT-E_POT2+VdcNd-DMFT.EDC_imp
         # Print out results
         main_iter.write( '%3d %3d %14.6f %10.6f %14.6f %10.6f' %(itt+1,it+1,ETOT,ETOT-ETOT_old,ETOT2,ETOT2-ETOT2_old) )
         main_iter.write('\n')
         main_iter.flush()
 
#      if TB.LHF=='.FALSE.': 
#         Tot_SigAvg=zeros(TB.ncor_orb)
#         for i,ats in enumerate(cor_at):
#            d_orb=TB.TB_orbs[ats[0]]
#            dir_name='imp.'+str(i)+'/'
#            SigAvg=Read_Prob.Compute_TrSigmaG(dir_name,U[i],TB.ncor_orb_list[i])
#            for at in ats: 
#               for orb in TB.TB_orbs[at]:
#                  idx=TB.idx[at][orb]
#                  Tot_SigAvg[idx] = SigAvg[d_orb.index(orb)] - DMFT.VDC[i]
#         fi=open('SigAvg.out','w')
#         for data in Tot_SigAvg: 
#            print >>fi, data,
#         print >>fi, ''
#         fi.close()
#         if p['nspin']==2:
#            fi=open('SigAvg_dn.out','w')
#            for data in DMFT.SigMdc_dn:
#               print >>fi, data,
#            print >>fi, ''
#            fi.close()
#      else: 
#         fi=open('SigAvg.out','w')
#         for data in DMFT.SigMdc:
#            print >>fi, data,
#         print >>fi, ''
#         fi.close()
#         if p['nspin']==2:
#            fi=open('SigAvg_dn.out','w')
#            for data in DMFT.SigMdc_dn:
#               print >>fi, data,
#            print >>fi, ''
#            fi.close()

      n_avg=int(p['Nit'])
      fiDMFT=open('ENERGY','r')
      Eline=fiDMFT.readlines()[-n_avg:]
      EAVG=0.; EAVG2=0.; 
      for i in range(n_avg):
         EAVG+=float(Eline[i].split()[2])
         EAVG2+=float(Eline[i].split()[4])
      fiDMFT.close()
      main_iter.write( '--------------------------\n' )
      main_iter.write( '%14.6f %14.6f %10.6f ' %(EAVG/n_avg, EAVG2/n_avg, CHGDIFF) )
      main_iter.write('\n')
      main_iter.write( '--------------------------\n' )
      main_iter.flush()
      E_iter.write( '%14.6f %14.6f' %(EAVG/n_avg, EAVG2/n_avg) )
      E_iter.write('\n')
      E_iter.flush()

   main_out.write( 'Caculation Ends'+now() )
   main_out.write('\n')
   main_out.flush()


