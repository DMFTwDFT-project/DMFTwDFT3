from scipy import *

def Compute_TrSigmaG(path,U,n_orb=5):
   fi=open(path+'actqmc.cix','r')
   for i in range(2): fi.readline()
   line=fi.readline()
   n_state=int(line.split()[1])
   n_bath=int(line.split()[2])
   if 2*n_orb!=n_bath: print "Something is wrong in actqmc.cix!!"; exit()
   #print n_state,n_bath  
   symidx=zeros(n_bath,dtype=int)
   e_orb=zeros(n_bath)
   fi.readline()
   for i in range(n_bath): symidx[i]=int(fi.readline().split()[2])
   fi.readline()
   e_orb_sym=map(float,fi.readline().split())
   for i in range(n_bath): e_orb[i]=e_orb_sym[symidx[i]]
   #print e_orb
   fi.readline()
   states=zeros((n_state,n_bath),dtype=int)
   for i in range(n_state): 
      line=fi.readline()
      #print line.split()[0], line.split()[-1] 
      occ=line.split()[-1] 
      for j,oc in enumerate(occ): 
         if oc=='0': pass
         elif oc=='u': states[i,j]=1
         elif oc=='d': states[i,j+n_orb]=1
         elif oc=='2': states[i,j]=1;states[i,j+n_orb]=1
         else: print "Something is wrong in actqmc.cix"; exit()
   #print states[0],states[1]
   fi.close()
   fi=open(path+'UC.dat','r')
   UC=[]
   for line in fi.readlines():
      UC.append(map(float,line.split()))
   if len(UC)!=2*n_orb: print "The size of UC is not consistent with orb"; exit()
   UC=array(UC)+U-diag(ones(2*n_orb)*U)
   #print UC
   Prob=zeros(n_state)
   fi=open(path+'Probability.dat','r')
   fi.readline()
   for i in range(n_state): Prob[i]=float(fi.readline().split()[2])
   fi.close()
   #print Prob
   Nd=0
   for i in range(n_state): Nd+=sum(states[i])*Prob[i]
   #print Nd
   mom=zeros(n_orb)
   for i in range(n_state): 
      for j in range(n_orb): mom[j]+=states[i][j]*Prob[i];mom[j]+=states[i][j+n_orb]*Prob[i]
   Eorb=0
   for i in range(n_state): Eorb+=e_orb.dot(states[i])*Prob[i]
   #print Eorb
   Epot=0
   for i in range(n_state): Epot+=0.5*states[i].dot(UC.dot(states[i]))*Prob[i]
   #print Epot+Eorb
   Sigma_avg=zeros(n_orb)
   for i in range(1,n_state):
      Epot2=states[i]*(UC.dot(states[i]))
      for j in range(n_orb): 
         Sigma_avg[j]+=(Epot2[j]+Epot2[j+n_orb])*Prob[i]
   return Sigma_avg/mom
   #print Sigma_avg/mom

if __name__=='__main__':

   Compute_TrSigmaG(5)
