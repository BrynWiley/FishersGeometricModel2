import csv
import os
from multiprocessing import Pool
import sys
#This file is used to enable multiple haploid SLiM simulations to be run in parallel,
#using Haploid_Finite_Population_Creation.slim.
#Designed to be run through the command line

#Define the l, lambda, library_id, and maximum number of parallel processes through the command line
l=int(float(sys.argv[1]))
lam=float(sys.argv[2])
library_id_start=int(float(sys.argv[3]))
max_processes=int(float(sys.argv[4]))

def runcommand(values):
  runNumber=values[0]
  library_id=values[1]
  k=values[2]
  mutationRate=values[3]
  command="/Linux/bin/slim -d runNumber={} -d library_id={} -d k={} -d mutationRate={} -d lambda={} -d genome_size={} -s {} Haploid_Finite_Population_Creation.slim".format(runNumber,library_id,k,mutationRate,lam,l,runNumber*library_id*k*l)
  os.system(command)
  
runNumber=1
library_id=library_id_start
command_list=[]

#Generate SLiM simulation commands for the range of parameter values required for this experiment
#mu: mutation rate
#k: the k value for fitness calculation

for replicate in range(25):
  for k in [2,6,10]:
    for mu in [0.1,0.01,0.001]:
      command=(runNumber,library_id,k,mu/l)
      command_list.append(command)
      runNumber=runNumber+1
  library_id=library_id+1

with Pool(processes=max_processes) as pool:
  pool.map(runcommand,command_list,1)



          
            
