import os, time, subprocess

if __name__=='__main__':
   trigger='THz'
   trig1='dx'
   num_at= int(input('how many atoms do you have?'))
   homer=range(1,num_at+1,1)
   flag = 0
   idof=1
   begin=10**10

#   holder=None

   if os.path.exists('OUTCAR') and os.path.exists('DYNMAT'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('I see you')
            begin=i
            flag+=1
#            kpt=1
#            if flag>0:
#               Dells=open('Ramode_%s.deig' %(str(flag)),'w')
#         if i>begin:
#            if trig1 in lines[i]:
            idof=1
            count=1
            for j in range(i+2,i+num_at+2):
               cx= float(lines[j].split()[3])
               cy= float(lines[j].split()[4])
               cz= float(lines[j].split()[5])
               cmd='awk \'{print $1, $2, $3*'+str(cx)+'}\' deltEps_'+str(idof)+'.deig > deltEps_1_'+str(idof)+'.deig'
               print os.popen(cmd).read()
               cmd='awk \'{print $1, $2, $3*'+str(cy)+'}\' deltEps_'+str(idof+1)+'.deig > deltEps_2_'+str(idof)+'.deig'
               print os.popen(cmd).read()
               cmd='awk \'{print $1, $2, $3*'+str(cz)+'}\' deltEps_'+str(idof+2)+'.deig > deltEps_3_'+str(idof)+'.deig'
               print os.popen(cmd).read()
               cmd='paste deltEps_{1,2,3}_'+str(idof)+'.deig | awk \'{ print $1, $2, $3+$6+$9 }\'>Ramode_'+str(flag)+'_'+str(count)+'.deig'
               print os.popen(cmd).read()
               idof+=3
               count+=1
#              cmd='paste Ramode_'+str(flag)+'_{
