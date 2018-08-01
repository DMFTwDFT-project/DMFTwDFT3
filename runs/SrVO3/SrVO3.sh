#!/bin/bash
#PBS -N SrVO3
#PBS -q alromero #comm_256g_mem 
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=16,pvmem=10gb
#PBS -m ae
#PBS -M ukh0001@mix.wvu.edu
#PBS -j oe
#PBS -e /users/ukh0001/projects/DFTDMFT/runs/SrVO3/OUTPUT.error
#PBS -o /users/ukh0001/projects/DFTDMFT/runs/SrVO3/OUTPUT.output

source ~/.bashrc
ulimit -s unlimited
cd  /users/ukh0001/projects/DFTDMFT/runs/SrVO3
time python RUNDMFT.py  
echo 'Done'

