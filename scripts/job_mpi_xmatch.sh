
#!/bin/bash
#PBS -V
#PBS -N xmatch
#PBS -k eo
#PBS -l nodes=2:ppn=40
#PBS -l walltime=100:00:00
#################################

### Switch to the working directory;
cd $PBS_O_WORKDIR

### Run:
module load openmpi-1.10.4
mpirun -np 80 python3 /home/mcavieres/halo_map/scripts/xmatch_MPI.py
echo "done"

