import os, time, subprocess

if __name__=='__main__':

   trigger='Degree'
   trig1='band No.'
   bands= int(input('how many bands do you have?'))
   flag = -1
   begin=10**6
   kpt=1

   if os.path.exists('OUTCAR'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('hello darkness my old friend')
            begin=i
            flag+=1
            kpt=1
            if flag>0:
               Dells=open('deltEps_%s.deig' %(str(flag)),'w')
         if i>begin:
            if trig1 in lines[i]:
               for j in range(i+1,i+1+bands):
                  bandno=int(lines[j].split()[0])
                  celta= float(lines[j].split()[1])
                  Dells.write("%d       %d      %.6f    \n" %(bandno,kpt,celta))
               kpt=kpt+1


   trigger='THz'
   trig1='dx'
   num_at= int(input('how many atoms do you have?'))
   flag= 0
   strang=''
   strangat=''

   for i in range(num_at):
      strang+=str(i+1)+','
   strang=strang[:-1]

#   strangtxt=open("strangalang",'w')
#   n=strangtxt.write(strang)
#   strangtxt.close()
   for i in range(num_at):
      strangat+='$'+str(3*(i+1))+'+'
   strangat=strangat[:-1]

#   strangattxt=open("strangatalang",'w')
#   n=strangattxt.write(strangat)
#   strangattxt.close()
#   print(strangat)

   if os.path.exists('OUTCAR') and os.path.exists('DYNMAT'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('I see you')
            flag+=1
            count=1
            idof=1
            for j in range(i+2,i+num_at+2):
               cx= float(lines[j].split()[3])
               cy= float(lines[j].split()[4])
               cz= float(lines[j].split()[5])
               cmd='awk \'{print $1, $2, $3*'+str(cx)+'}\' deltEps_'+str(idof)+'.deig > deltEps_1_'+str(flag)+'_'+str(count)+'.deig'
               print os.popen(cmd).read()
               cmd='awk \'{print $1, $2, $3*'+str(cy)+'}\' deltEps_'+str(idof+1)+'.deig > deltEps_2_'+str(flag)+'_'+str(count)+'.deig'
               print os.popen(cmd).read()
               cmd='awk \'{print $1, $2, $3*'+str(cz)+'}\' deltEps_'+str(idof+2)+'.deig > deltEps_3_'+str(flag)+'_'+str(count)+'.deig'
               print os.popen(cmd).read()
               cmd='paste deltEps_{1,2,3}_'+str(flag)+'_'+str(count)+'.deig | awk \'{ print $1, $2,$3+$6+$9 }\'>Ramode_'+str(flag)+'_'+str(count)+'.deig'
               print os.popen(cmd).read()
               idof+=3
               count+=1
            cmd='paste Ramode_'+str(flag)+'_{'+strang+'}.deig | awk \'{print $1,$2,'+strangat+'}\'>Ramode_'+str(flag)+'.deig'
            print os.popen(cmd).read()

      cmd='rm Ramode_*_*.deig | rm deltEps_*.deig'
      print os.popen(cmd).read()
                                                                                                                                                                                      81,0-1        Bo
