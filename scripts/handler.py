import os, time, subprocess

if __name__=='__main__':
   sez=str(raw_input("What dimensions do you want your supercell?\n Enter 3 integers with one whitespace between them and backslash double-quotes (e.g \\\") around all three: "))
   cmd="phonopy -d --dim="+sez

   print os.popen(cmd).read()

   for ct in range(3):
      cmd="mv POSCAR POSCAR_nawt"
      print os.popen(cmd).read()
      cmd="mkdir VASP_FD_"+str(ct+1)
      print os.popen(cmd).read()
      cmd = "cp ./* VASP_FD_"+str(ct+1)
      print os.popen(cmd).read()
      cmd="cd VASP_FD_"+str(ct+1)
      print os.popen(cmd).read()
      cmd="mv POSCAR-00"+str(ct+1)
      print os.popen(cmd).read()
      cmd="cd .."

