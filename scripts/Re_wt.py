import os, time, subprocess

if __name__=='__main__':
   trigger=str(raw_input("Which deg. of fredom do you want from OUTCAR?\n (Please mind your whitespace and if mode is imaginary):"))

   count_me = int(raw_input("How many atoms are there:"))

   trig1='band No.'
   bands= 60
   dof = 1

   xwt=0.
   ywt=0.
   zwt=0.
   

   begin=10**10
   kpt=1
#   holder=None
   
   if os.path.exists('OUTCAR'):
      lines=open('OUTCAR','r').readlines()[1:]
      for i, line in enumerate(lines):
         if trigger in lines[i]:
            print('hello darkness my old friend')
            begin=i
            for el in range(i+2,i+2+count_me):
               xwt=float(lines[el].split()[3])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(xwt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               print os.popen(cmd).read()
               print(xwt) 
               dof+=1
               ywt=float(lines[el].split()[4])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(ywt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               print os.popen(cmd).read()
               dof+=1
               print(ywt)
               zwt=float(lines[el].split()[5])
               tmp_nome=str('deltEps_%s.deig' %str(dof))
               cmd="awk '{print $1, $2, $3*(%f)}' " %(zwt)+tmp_nome+" > deltEps_%s_primo.deig" %str(dof)
               print os.popen(cmd).read()
               dof+=1
               print(zwt)
#  
#            if flag>0: 
#               Dells=open('deltEps_%s.deig' %(str(flag)),'w')
#         if i>begin:
#            if trig1 in lines[i]:
#               for j in range(i+1,i+1+bands):
##                 holder=str( lines[j])
#                  bandno=int(lines[j].split()[0])
#                  celta= float(lines[j].split()[1])
#                 
#                  Dells.write("%d        %d        %.6f\n" %(bandno,kpt,celta))
#               kpt+=1
