#!/bin/bash -ex
#SBATCH -p hsw_p100
#SBATCH -N 1
##SBATCH -n 1
##SBATCH --qos=short
#SBATCH -t 06:00:00
#SBATCH -J MPAS
#SBATCH -o MPAS.out

# Load saved modules
module load null
module load PrgEnv/PGI+OpenMPI/2017-05-24
module load openmpi

export PATH=/cm/extra/apps/PGI/17.5/linux86-64/17.5/bin:/cm/extra/apps/PGI/17.5/linux86-64/2017/mpi/openmpi-1.10.2/bin:$PATH
export LD_LIBRARY_PATH=/cm/extra/apps/PGI/17.5/linux86-64/17.5/lib:/cm/extra/apps/PGI/17.5/linux86-64/2017/mpi/openmpi-1.10.2/lib:$LD_LIBRARY_PATH

make
