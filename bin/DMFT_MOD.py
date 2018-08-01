#!/usr/bin/env python

from scipy import *
import os,sys,copy
import Fileio
import scipy.interpolate

class DMFT_class:
   """
   A general class to manipulate DMFT modules
   """
   def __init__(self,p,pC,TB):

      cor_at=p['cor_at']
      cor_orb=p['cor_orb']
      self.cor_at=cor_at
      self.cor_orb=cor_orb
      self.mu=0.0
      self.EPOT=0.0
      self.EPOT2=0.0
      self.EDC=0.0
      self.EKIN=0.0
      self.EKIN0=0.0
      self.Natom=zeros(len(self.cor_at),dtype=int)
      self.Nd_latt=zeros(len(self.cor_at),dtype=float)
      self.VDC=zeros(len(self.cor_at),dtype=float)
      self.T=1.0/pC['beta'][0]

      self.N_latt=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      self.MOM=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      #for i in range(len(cor_at)): self.VDC[i]=(p['J'][i]-p['U'][i])/2.

      self.nom=p['nom']
      noms=p['noms']
      nomlog=p['nomlog']
     
      self.EqLogMesh(p['noms'],p['nomlog'],p['nom'],self.T)

      for i,ats in enumerate(cor_at):
         self.Natom[i]=len(cor_at[i])


   def Mix_Sig(self,nspin,mix_sig):
      self.Sig=array(self.Sig_old)+mix_sig*(array(self.Sig)-array(self.Sig_old))
      self.SigMdc=array(self.SigMdc_old)+mix_sig*(array(self.SigMdc)-array(self.SigMdc_old))
      if nspin==2:
         self.Sig_dn=array(self.Sig_dn_old)+mix_sig*(array(self.Sig_dn)-array(self.Sig_dn_old))
         self.SigMdc_dn=array(self.SigMdc_dn_old)+mix_sig*(array(self.SigMdc_dn)-array(self.SigMdc_dn_old))
      
   def Read_Sig(self,TB,nspin):
      self.Sigoo=[];SigMoo=[]
      self.Nd_imp=zeros(len(self.cor_at),dtype=float)
      self.N_imp=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      self.MOM_imp=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      self.Eimp=[]
      for i,ats in enumerate(self.cor_at):
         d_orb=TB.TB_orbs[ats[0]]
         fileSig='Sig'+str(i+1)+'.out'
         if (os.path.exists(fileSig)): # If output file exists, start from previous iteration
            (ommesh_long,Sig_file,TrS,Epot,nf_q,mom) = Fileio.Read_complex_Data(fileSig)
            if len(Sig_file)!=nspin*len(self.cor_orb[i]): print "The number of correated orbital is not same as Sig file column"; exit()
            if len(mom)!=nspin*len(self.cor_orb[i]): print "The number of correated orbital is not same as mom list in Sig file"; exit()
            if self.nom>len(ommesh_long): print "nom should be decreased!"; exit()
            self.Nd_imp[i]=nf_q
            for j,orbs in enumerate(self.cor_orb[i]):
               for orb in orbs:
                  self.N_imp[i,d_orb.index(orb)]=mom[j]/len(orbs)
                  #self.N_imp[TB.idx[at][orb]]=mom[j]/len(orbs)
                  if nspin==2: 
                     self.N_imp[i,d_orb.index(orb)]+=mom[j+len(self.cor_orb[i])]/len(orbs)
                     self.MOM_imp[i,d_orb.index(orb)]=(mom[j]-mom[j+len(self.cor_orb[i])])/len(orbs)
            #self.Eimp.append(TrS-0.5*sum([mom[j]*Sig_file[j][-1].real for j in range(len(mom))])); #nf_qmc.append(list(mom))
            self.Eimp.append(TrS); #nf_qmc.append(list(mom))
            self.Sigoo.append(array(Sig_file[:,-1].real))
            for iom in range(len(ommesh_long)):
               Sig_file[:,iom]-=self.Sigoo[i]
               for ii in range(nspin*len(self.cor_orb[i])):
                  if Sig_file[ii,iom].imag >0: Sig_file[ii,iom]=Sig_file[ii,iom].real+0.0j
            SigMoo.append(Interpolate(ommesh_long,Sig_file,self.ommesh,1))
         else:
            self.Eimp.append(0.0)
            self.Sigoo.append([])
            SigMoo.append([])
            for j in range(nspin*len(self.cor_orb[i])):
               if j<len(self.cor_orb[i]):
                  self.Sigoo[i].append(-0.1)
               else: self.Sigoo[i].append(0.1) #Break symmetry
               SigMoo[i].append(zeros(len(self.ommesh)))

      self.Sig=zeros((len(TB.cor_idx),len(self.ommesh)),dtype=complex)
      for i,ats in enumerate(self.cor_at):
         for at in ats:
            for j,orbs in enumerate(self.cor_orb[i]):
               for orb in orbs:
                  idx=TB.idx[at][orb]
                  self.Sig[idx,:] = copy.deepcopy(SigMoo[i][j])
      if (os.path.exists('SigMoo.out')):
         omtemp,self.Sig_old=Fileio.Read_complex_multilines('SigMoo.out')
      else: 
         self.Sig_old=copy.deepcopy(self.Sig)

      if nspin==2:
         self.Sig_dn=zeros((len(TB.cor_idx),len(self.ommesh)),dtype=complex)
         for i,ats in enumerate(self.cor_at):
            for at in ats:
               for j,orbs in enumerate(self.cor_orb[i]):
                  for orb in orbs:
                     idx=TB.idx[at][orb]
                     self.Sig_dn[idx,:] = copy.deepcopy(SigMoo[i][j+len(self.cor_orb[i])])
         if (os.path.exists('SigMoo_dn.out')):
            omtemp,self.Sig_dn_old=Fileio.Read_complex_multilines('SigMoo_dn.out')
         else: 
            self.Sig_dn_old=copy.deepcopy(self.Sig_dn)

   def Cmp_Sig_highT(self,T_high,noms_high,nspin):
      self.ommesh_highT=zeros(noms_high)
      for i in range(noms_high): self.ommesh_highT[i]=pi*T_high*(2*i+1)
      self.Sig_highT=Interpolate(self.ommesh,self.Sig,self.ommesh_highT,1)
      if nspin==2: self.Sig_dn_highT=Interpolate(self.ommesh,self.Sig_dn,self.ommesh_highT,1)
         
      
   def Compute_HF(self,Nd_qmc,p,TB):
      """Compute Energy and Sigma"""
      nspin=p['nspin'];U=p['U'];J=p['J'];dc_type=p['dc_type'];Uprime=p['Uprime']
      self.VHF=zeros((len(self.cor_at),2*TB.max_cor_orb),dtype=float)
      self.VHF_latt=zeros((len(self.cor_at),2*TB.max_cor_orb),dtype=float)
      self.VHF_imp=zeros((len(self.cor_at),2*TB.max_cor_orb),dtype=float)
      self.VDC=zeros(len(self.cor_at),dtype=float)
      self.sig_st=[]#zeros(len(cor_at),dtype=float)
      self.EHF=0.0#zeros(len(self.cor_at),dtype=float)
      self.EHF_latt=0.0#zeros(len(self.cor_at),dtype=float)
      self.EHF_imp=0.0#zeros(len(self.cor_at),dtype=float)
      self.EDC=0.0;self.EDC_imp=0.0
      #self.EHF_cor=zeros(len(self.cor_at),dtype=float)
      self.SigMdc=zeros(TB.ncor_orb,dtype=float)
      if nspin==2: self.SigMdc_dn=zeros(TB.ncor_orb,dtype=float)
      for i,ats in enumerate(self.cor_at):
         d_orb=TB.TB_orbs[ats[0]]
         fi=open('UC'+str(i+1)+'.dat','r')
         UC=[]
         for line in fi.readlines():
            UC.append(map(float,line.split()))
         if len(UC)!=2*len(d_orb): print "The size of UC is not consistent with orb"; exit()
         UC=array(UC)+U[i]-diag(ones(2*len(d_orb))*U[i])
         if Nd_qmc==0:
            OCC=array(list((self.N_latt[i][:len(d_orb)]+self.MOM[i][:len(d_orb)])/2)+list((self.N_latt[i][:len(d_orb)]-self.MOM[i][:len(d_orb)])/2))
         else:
            OCC=array(list((self.N_imp[i][:len(d_orb)]+self.MOM_imp[i][:len(d_orb)])/2)+list((self.N_imp[i][:len(d_orb)]-self.MOM_imp[i][:len(d_orb)])/2))
         OCC_latt=array(list((self.N_latt[i][:len(d_orb)]+self.MOM[i][:len(d_orb)])/2)+list((self.N_latt[i][:len(d_orb)]-self.MOM[i][:len(d_orb)])/2))
         OCC_imp=array(list((self.N_imp[i][:len(d_orb)]+self.MOM_imp[i][:len(d_orb)])/2)+list((self.N_imp[i][:len(d_orb)]-self.MOM_imp[i][:len(d_orb)])/2))
         self.VHF[i,:2*len(d_orb)]=dot(UC,OCC)
         self.VHF_latt[i,:2*len(d_orb)]=dot(UC,OCC_latt)
         self.VHF_imp[i,:2*len(d_orb)]=dot(UC,OCC_imp)

         ###### Compute VDC #######
         if dc_type==1:
            if Nd_qmc==0:
               self.VDC[i]=Uprime[i]*(self.Nd_latt[i]-0.5)-J[i]/2*(self.Nd_latt[i]-1)
            else:
               self.VDC[i]=Uprime[i]*(self.Nd_imp[i]-0.5)-J[i]/2*(self.Nd_imp[i]-1)

         elif dc_type==2:
            if Nd_qmc==0:
               self.VDC[i]=(Uprime[i]-2*J[i])*(0.9*self.Nd_latt[i])-J[i]/2*(2*self.Nd_latt[i]/5)
            else:
               self.VDC[i]=(Uprime[i]-2*J[i])*(0.9*self.Nd_imp[i])-J[i]/2*(2*self.Nd_imp[i]/5)
         else: print "dc type is wrong!"; exit()

         self.sig_st.append([])
         for j,orbs in enumerate(self.cor_orb[i]):
            self.sig_st[i].append(0.0)
            for orb in orbs:
               idx=d_orb.index(orb)
               #for orb2 in orbs:
               #   idx2=d_orb.index(orb2)
               #   ########### This only works for one cor_orb #############
               #   self.EHF_cor[i]+=0.5*(OCC[idx]*UC[idx,idx2]*OCC[idx2]+OCC[idx]*UC[idx,idx2+5]*OCC[idx2+5]+OCC[idx+5]*UC[idx+5,idx2]*OCC[idx2]+OCC[idx+5]*UC[idx+5,idx2+5]*OCC[idx2+5])
               for orb2 in TB.TB_orbs[ats[0]]:
                  if TB.cor_idx[TB.idx[ats[0]][orb2]]==1:
                     idx2=d_orb.index(orb2)
                     self.sig_st[i][j]+=UC[idx,idx2]*OCC[idx2]+UC[idx,idx2+len(d_orb)]*OCC[idx2+len(d_orb)]
            self.sig_st[i][j]/=len(orbs)
            for orb in orbs:
               idx=d_orb.index(orb)
               if nspin==1 and Nd_qmc>1: self.VHF[i,idx]=self.Sigoo[i][j]; self.VHF[i,idx+len(d_orb)]=self.Sigoo[i][j]
               if nspin==2 and Nd_qmc>1: self.VHF[i,idx]=self.Sigoo[i][j]; self.VHF[i,idx+len(d_orb)]=self.Sigoo[i][j+len(self.cor_orb[i])]
            self.sig_st[i][j]-=self.VDC[i]
         self.EHF+=len(ats)*0.5*dot(OCC,self.VHF[i][:2*len(d_orb)])
         self.EHF_imp+=len(ats)*0.5*dot(OCC_imp,self.VHF_imp[i][:2*len(d_orb)])
         self.EHF_latt+=len(ats)*0.5*dot(OCC_latt,self.VHF_latt[i][:2*len(d_orb)])

         if dc_type==1:
            self.EDC+=len(ats)*(Uprime[i]*self.Nd_latt[i]*(self.Nd_latt[i]-1)/2.0-J[i]*self.Nd_latt[i]*(self.Nd_latt[i]-2)/4.0)
            self.EDC_imp+=len(ats)*(Uprime[i]*self.Nd_imp[i]*(self.Nd_imp[i]-1)/2.0-J[i]*self.Nd_imp[i]*(self.Nd_imp[i]-2)/4.0)
         elif dc_type==2:
            self.EDC+=self.Natom[i]*(Uprime[i]*self.Nd_latt[i]*(0.9*self.Nd_latt[i])/2-5*J[i]*self.Nd_latt[i]/2*(2*self.Nd_latt[i]/5))
            self.EDC_imp+=self.Natom[i]*(Uprime[i]*self.Nd_imp[i]*(0.9*self.Nd_imp[i])/2-5*J[i]*self.Nd_imp[i]/2*(2*self.Nd_imp[i]/5))
         else: print "This dc type is not supported!"; exit()

         for at in ats:
            for orb in TB.TB_orbs[at]:
               idx=TB.idx[at][orb]
               self.SigMdc[idx] = self.VHF[i,d_orb.index(orb)]-self.VDC[i]
         if nspin==2: 
            for at in ats:
               for orb in TB.TB_orbs[at]:
                  idx=TB.idx[at][orb]
                  self.SigMdc_dn[idx] = self.VHF[i,d_orb.index(orb)+len(d_orb)]-self.VDC[i]
      if (os.path.exists('SigMdc.out')):
         self.SigMdc_old=Fileio.Read_float('SigMdc.out')
      else: 
         self.SigMdc_old=copy.deepcopy(self.SigMdc)
      if nspin==2:
         if (os.path.exists('SigMdc_dn.out')):
            self.SigMdc_dn_old=Fileio.Read_float('SigMdc_dn.out')
         else:
            self.SigMdc_dn_old=copy.deepcopy(self.SigMdc_dn)


   def Compute_Delta(self,T,nspin,cor_at,cor_orb,TB,nom,delta=0.0):
   #####  Store local Green function and local Self energy with equidistant mesh as a list type ##########
      ommesh,tot_GLOC=Fileio.Read_complex_multilines('G_loc.out')    
      ommesh,tot_SIGLOC=Fileio.Read_complex_multilines('SigMoo.out')    
      SIGMDC=loadtxt('SigMdc.out')
      if nspin==2:
         ommesh,tot_GLOC_dn=Fileio.Read_complex_multilines('G_loc_dn.out')
         ommesh,tot_SIGLOC_dn=Fileio.Read_complex_multilines('SigMoo_dn.out')
         SIGMDC_dn=loadtxt('SigMdc_dn.out')
      DMFT_mu=loadtxt('DMFT_mu.out')
      tot_ed=loadtxt('Ed.out')
#      tot_sig_st=loadtxt('Sig_st.out')
      ed=[]; 
      GLOC=[];  SIGLOC=[]
      for i,ats in enumerate(cor_at):
         GLOC.append([])
         SIGLOC.append([])
         ed.append([])
         #sig_st.append([])
         for j,orbs in enumerate(cor_orb[i]):
            Gf_avg=zeros(len(ommesh),dtype=complex)
            Sig_avg=zeros(len(ommesh),dtype=complex)
            ed_avg=0.0
            #sigst_avg=0.0
            for at in ats:
               for orb in orbs:
                  idx=TB.idx[at][orb]
                  Gf_avg+=tot_GLOC[idx];Sig_avg+=tot_SIGLOC[idx]+SIGMDC[idx]
                  ed_avg+=tot_ed[idx]
                  #sigst_avg+=DMFT.Sig_st[i][j]#tot_sig_st[idx]
            Gf_avg/=len(ats)*len(orbs)
            Sig_avg/=len(ats)*len(orbs)#;Sig_avg-=sig_st[i]
            ed_avg/=len(ats)*len(orbs)
            #sigst_avg/=len(ats)*len(orbs)
            GLOC[i].append(list(Gf_avg))
            SIGLOC[i].append(list(Sig_avg))
            ed[i].append(ed_avg)
            #sig_st[i].append(sigst_avg)
      if nspin==2:
         GLOC_dn=[];  SIGLOC_dn=[]
         for i,ats in enumerate(cor_at):
            GLOC_dn.append([])
            SIGLOC_dn.append([])
            for j,orbs in enumerate(cor_orb[i]):
               Gf_avg=zeros(len(ommesh),dtype=complex)
               Sig_avg=zeros(len(ommesh),dtype=complex)
               for at in ats:
                  for orb in orbs:
                     idx=TB.idx[at][orb]
                     Gf_avg+=tot_GLOC_dn[idx];Sig_avg+=tot_SIGLOC_dn[idx]+SIGMDC_dn[idx]
               Gf_avg/=len(ats)*len(orbs)
               Sig_avg/=len(ats)*len(orbs)#;Sig_avg-=sig_st[i]
               GLOC_dn[i].append(list(Gf_avg))
               SIGLOC_dn[i].append(list(Sig_avg))
      for i in range(len(GLOC)):
         if len(GLOC[i])>0:
            Delta_s=zeros((nspin*len(GLOC[i]),len(ommesh)),dtype=complex)
            for j in range(len(GLOC[i])):
               Delta_s[j,:]=1j*ommesh+DMFT_mu-ed[i][j]-array(SIGLOC[i][j])+1j*delta-1.0/array(GLOC[i][j])
            if nspin==2:
               for j in range(len(GLOC[i])):
                  Delta_s[j+len(GLOC[i]),:]=1j*ommesh+DMFT_mu-ed[i][j]-array(SIGLOC_dn[i][j])+1j*delta-1.0/array(GLOC_dn[i][j])
         ######  Interpolate Delta ####
         ommesh_new=pi*T*(2*arange(nom)+1)
   
         Delta=Interpolate(ommesh,Delta_s,ommesh_new,1)
         Fileio.Print_complex_multilines(Delta,ommesh_new,'Delta'+str(i+1)+'.inp')
      return DMFT_mu,ed#,sig_st


   def mu_bcast(self,comm):
      self.mu = comm.bcast(self.mu, root=0)

   def Gather_Ksum_HF(self,comm,TOT_NKPTS,NWANN):
      rank=comm.Get_rank()
      size=comm.Get_size()
      self.EKIN=comm.gather(self.EKIN, root=0)
      if rank==0: self.EKIN=sum(self.EKIN)/TOT_NKPTS*2
      if rank==0: dm=zeros((size,NWANN),dtype=float)
      else: dm=None
      comm.Gather(self.dm, dm, root=0)
      if rank==0:
         self.dm=zeros(NWANN,dtype=float)
         for i in range(NWANN):
            self.dm[i]=sum(dm[:,i])/TOT_NKPTS*2


   def Gather_Ksum_DMFT(self,comm,TOT_NKPTS,NWANN):
      rank=comm.Get_rank()
      size=comm.Get_size()
      self.EKIN=comm.gather(self.EKIN, root=0)
      if rank==0: self.EKIN=sum(self.EKIN)/TOT_NKPTS*2
      if rank==0: dm=zeros((size,NWANN),dtype=float)
      else: dm=None
      comm.Gather(self.dm, dm, root=0)
      if rank==0:
         self.dm=zeros(NWANN,dtype=float)
         for i in range(NWANN):
            self.dm[i]=sum(dm[:,i])/TOT_NKPTS*2
      if rank==0: Gloc=zeros((size,self.ncor_orb,len(self.ommesh)),dtype=complex)
      else: Gloc=None
      comm.Gather(self.Gloc, Gloc, root=0)
      if rank==0:
         self.Gloc=zeros((self.ncor_orb,len(self.ommesh)),dtype=complex)
         for i in range(self.ncor_orb):
            for iom in range(len(self.ommesh)):
               self.Gloc[i,iom]=sum(Gloc[:,i,iom])/TOT_NKPTS

   def Gather_Ksum_DMFT2(self,comm,NFDIR,TOT_NKPTS,NWANN):
      rank=comm.Get_rank()
      size=comm.Get_size()
      if rank==0: ddm=zeros((size,NFDIR,NWANN),dtype=float)
      else: ddm=None
      comm.Gather(self.ddm, ddm, root=0)
      if rank==0:
         self.ddm=zeros((NFDIR,NWANN),dtype=float)
         for i in range(size):
            self.ddm[:,:]+=ddm[i,:,:]
         self.ddm=self.ddm/TOT_NKPTS*2
      if rank==0: dGloc=zeros((size,NFDIR,self.ncor_orb,len(self.ommesh)),dtype=complex)
      else: dGloc=None
      comm.Gather(self.dGloc, dGloc, root=0)
      if rank==0:
         self.dGloc=zeros((NFDIR,self.ncor_orb,len(self.ommesh)),dtype=complex)
         for i in range(size):
            self.dGloc[:,:,:]+=dGloc[i,:,:,:]
         self.dGloc=self.dGloc/TOT_NKPTS

   def EqLogMesh(self,noms,nomlog,nom,T):
      """ This function computes the linear mesh for samll omega and the log mesh for large omega
      """
      self.ommesh=zeros(noms+nomlog,dtype=float)
      for i in range(noms): self.ommesh[i]=pi*T*(2*i+1)
      logi=log(pi*T*(2*noms+1));logf=log(pi*T*(2*(nom-1)+1))
      for i in range(noms,noms+nomlog): self.ommesh[i]=exp(logi+(logf-logi)*(i-noms)/(nomlog-1))

   def Update_Nlatt(self,TB,nspin=2):
      self.N_latt=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      self.MOM=zeros((len(self.cor_at),TB.max_cor_orb),dtype=float)
      self.Nd_latt=zeros(len(self.cor_at),dtype=float)
      line=open('DMINFO','r').readlines()[-1].split()
      for i,ats in enumerate(self.cor_at):
         d_orb=TB.TB_orbs[ats[0]]
         for at in ats:
            for j,orb in enumerate(d_orb):
               idx=TB.idx[at][orb]
               self.N_latt[i,j]+=float(line[idx])
               if nspin==2: self.MOM[i,j]+=float(line[idx+TB.ncor_orb])
         self.Nd_latt[i]=sum(self.N_latt[i])
         self.N_latt[i,:]/=len(ats)
         if nspin==2: self.MOM[i,:]/=len(ats)
         self.Nd_latt[i]/=len(ats)
   

   def Compute_VDC_latt(self,Uprime,J,dc_type):
      self.VDC_old=copy.deepcopy(self.VDC)
      self.VDC=zeros(len(self.Nd_latt),dtype=float)
      for i in range(len(self.Nd_latt)):
         if dc_type==1:
            #self.VDC[i]=(Uprime[i]-2*J[i])*(self.Nd_latt[i]-0.5)-J[i]/2*(self.Nd_latt[i]-3)
            self.VDC[i]=Uprime[i]*(self.Nd_latt[i]-0.5)-J[i]/2*(self.Nd_latt[i]-1)
            #self.VDC[i]=Uprime[i]*(self.N_imp[i]-0.5)-J[i]/2*(self.N_imp[i]-1)
         elif dc_type==2:
            self.VDC[i]=(Uprime[i]-2*J[i])*(0.9*self.Nd_latt[i])-J[i]/2*(2*self.Nd_latt[i]/5)
         else: print "dc type is wrong!"; exit()

   def Mix_DC(self,Uprime,J,self_dc,mix_dc,Nd_f):
      for i in range(len(self.Nd_latt)): 
         if self_dc==True:
            self.VDC[i]=self.VDC_old[i]+mix_dc*(self.VDC[i]-self.VDC_old[i])
         else:
            VDC_f=(Uprime[i]-2*J[i])*(Nd_f[i]-0.5)-J[i]/2*(Nd_f[i]-3)
            self.VDC[i]=self.VDC_old[i]+mix_dc*(VDC_f-self.VDC[i])

    
   def dm_update(self,cor_at):
      fi=open('dm.out','r')
      self.mu=float(fi.readline().split()[1])
      dc_array=fi.readline().split()
      if len(dc_array)==len(cor_at)+1:
         for i in range(len(cor_at)):
            self.VDC[i]=float(dc_array[i+1])
      else:
         for i in range(len(cor_at)):
            self.VDC[i]=float(dc_array[1])
      for i in range(3): fi.readline()
      self.dm=array([float(dmi) for dmi in fi.readline().split()])
      if len(self.VDC)!=len(cor_at): print "Specify VDC in dm.out as many as correlated atom"; exit()
      fi.close()
  

   def Mix_SigMdc(self,mix_sig):
      #self.SigMdc=self.SigMdc_old+mix_sig*(self.SigMdc-self.SigMdc_old)
      self.SigMdc=array(self.SigMdc_old)+mix_sig*(array(self.SigMdc)-array(self.SigMdc_old))


   def Compute_FORCE(self,cor_at,cor_orb,TB,nom,NFDIR):
      """ This function computes the FORCE of DMFT using E_pot=1/2*Tr[G_loc*Sig_loc]
          The static part Estat is HF part (omega->infty) and dynamtical part is DMFT part
      """
      self.FORCE=zeros(NFDIR)
      for ifd,ifdir in enumerate(range(NFDIR)):
         self.dGLOC=[]
         for i,ats in enumerate(cor_at):
            self.dGLOC.append([])
            for j,orbs in enumerate(cor_orb[i]):
               dGf_avg=zeros(len(self.ommesh),dtype=complex)
               for at in ats:
                  for orb in orbs:
                     idx=TB.idx[at][orb]
                     dGf_avg+=self.dGloc[ifd][idx]
               dGf_avg/=len(ats)*len(orbs)
               self.dGLOC[i].append(list(dGf_avg))
         ######  Interpolate ####
         ommesh_new=pi*self.T*(2*arange(nom)+1)
         dGf_loc=[]; Sig_loc=[]
         for i in range(len(self.dGLOC)):
            dGf_loc.append([]); Sig_loc.append([])
            Gf_intpol=Interpolate(self.ommesh,self.dGLOC[i],ommesh_new,1)
            Sig_intpol=Interpolate(self.ommesh,self.SIGLOC[i],ommesh_new,1)
            for j in range(len(self.dGLOC[i])):
               dGf_loc[i].append(list(Gf_intpol[j]))
               Sig_loc[i].append(list(Sig_intpol[j]))



         d_orb=['d_z2','d_x2y2','d_xz','d_yz','d_xy']
         loc_idx=0
         Fstat=0.0
         Fdyna=0.0
         for i,ats in enumerate(cor_at):
            for at in ats:
               for orb in d_orb:
                  idx=TB.idx[at][orb]
                  Fstat+=self.SigMdc[idx,-1].real*self.ddm[ifd,idx]
               ####  DMFT part ################
               # factor 4 = negative omega, spin
            for j,orbs in enumerate(cor_orb[i]):
               Sig0=Sig_loc[i][j][-1].real;Sig1=Sig_loc[i][j][-1].imag*ommesh_new[-1];Gf1=-dGf_loc[i][j][-1].imag*ommesh_new[-1];
               Fdyna+=len(ats)*len(orbs)*((sum(array(dGf_loc[i][j])*(array(Sig_loc[i][j])-Sig0)).real-Gf1*Sig1*sum(1/ommesh_new**2))*self.T*4+Gf1*Sig1/2/self.T)
         self.FORCE[ifd]=Fstat+Fdyna
         #print Fstat,Fdyna


   def Compute_EPOT(self,cor_at,cor_orb,TB,nom):
      """ This function computes the POT energy of DMFT using E_pot=1/2*Tr[G_loc*Sig_loc]
          The static part Estat is HF part (omega->infty) and dynamtical part is DMFT part
      """
      ######  Interpolate ####
      ommesh_new=pi*self.T*(2*arange(nom)+1)
      Gf_loc=[]; Sig_loc=[]
      for i in range(len(self.GLOC)):
         Gf_loc.append([]); Sig_loc.append([])
         Gf_intpol=Interpolate(self.ommesh,self.GLOC[i],ommesh_new,1)
         Sig_intpol=Interpolate(self.ommesh,self.SIGLOC[i],ommesh_new,1)
         for j in range(len(self.GLOC[i])):
            Gf_loc[i].append(list(Gf_intpol[j]))
            Sig_loc[i].append(list(Sig_intpol[j]))
         #for j in range(len(self.GLOC[i])):
         #   Gf_intpol=zeros(nom,dtype=complex)
         #   Sig_intpol=zeros(nom,dtype=complex)
         #   GfSpline = scipy.interpolate.splrep(ommesh, array(Gf_loc_small[i][j]).real, k=3, s=0)
         #   Gf_intpol += scipy.interpolate.splev(ommesh_new, GfSpline)
         #   GfSpline = scipy.interpolate.splrep(ommesh, array(Gf_loc_small[i][j]).imag, k=3, s=0)
         #   Gf_intpol += 1j*scipy.interpolate.splev(ommesh_new, GfSpline)
         #   Gf_loc[i].append(list(Gf_intpol))
         #   SigSpline = scipy.interpolate.splrep(ommesh, array(Sig_loc_small[i][j]).real, k=3, s=0)
         #   Sig_intpol += scipy.interpolate.splev(ommesh_new, SigSpline)
         #   SigSpline = scipy.interpolate.splrep(ommesh, array(Sig_loc_small[i][j]).imag, k=3, s=0)
         #   Sig_intpol += 1j*scipy.interpolate.splev(ommesh_new, SigSpline)
         #   Sig_loc[i].append(list(Sig_intpol))



      self.EPOT=0.0;Estat=0.0;Estat2=0.0;Edyna=0.0;self.EPOT2=0.0;Edyna2=0.0

      d_orb=['d_z2','d_x2y2','d_xz','d_yz','d_xy']
      loc_idx=0
      for i,ats in enumerate(cor_at):
         Estat+=len(ats)*self.EHF[i]
         Estat2+=len(ats)*self.EHF_cor[i]
         #for at in ats:
            #for orb in d_orb:
               #idx=TB.idx[at][orb]
               #Estat+=0.5*(self.SigMdc[idx,-1].real+self.VDC[i])*self.dm[idx]
            ####  DMFT part ################
            # factor 4 = negative omega, spin
         for j,orbs in enumerate(cor_orb[i]):
            Sig0=Sig_loc[i][j][-1].real;Sig1=Sig_loc[i][j][-1].imag*ommesh_new[-1]
            Edyna+=len(ats)*len(orbs)*0.5*((sum(array(Gf_loc[i][j])*(array(Sig_loc[i][j])-Sig0)).real-Sig1*sum(1/ommesh_new**2))*self.T*4+Sig1/2/self.T)
         if len(cor_orb[i])>0: Edyna2+=len(ats)*self.Eimp[i]
      self.EPOT=Estat+Edyna
      self.EPOT2=Edyna2
      #self.EPOT2=Estat-Estat2+Edyna2


   def Print_Gloc(self,print_at,TB):
      for at in print_at:
         Fileio.Print_complex_multilines(array([self.Gloc[i] for i in range(TB.idx[at][TB.TB_orbs[at][0]],TB.idx[at][TB.TB_orbs[at][-1]]+1)]),self.ommesh,'G_loc_'+at+'.out')


def Interpolate(ommesh_input,data_input,ommesh_output,sorder):  

   data_output=zeros((shape(data_input)[0],len(ommesh_output)),dtype=complex)
   data_input=array(data_input)
   for i in range(len(data_input)):
      SSpline = scipy.interpolate.splrep(ommesh_input, data_input[i].real, k=sorder, s=0)
      data_output[i,:] += scipy.interpolate.splev(ommesh_output, SSpline)
      SSpline = scipy.interpolate.splrep(ommesh_input, data_input[i].imag, k=sorder, s=0)
      data_output[i,:] += 1j*scipy.interpolate.splev(ommesh_output, SSpline)
   return data_output
