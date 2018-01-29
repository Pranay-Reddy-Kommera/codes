#!/bin/tcsh
#
# LSF batch script to run an MPI application
#
#BSUB -P NCIS0002           # project code
#BSUB -W 01:00               # wall-clock time (hrs:mins)
#BSUB -n 1                  # number of tasks in job         
#BSUB -R "span[ptile=1]"    # run 16 MPI tasks per node
#BSUB -J dg2d                # job name
#BSUB -o ACCBaseCode1.out  # output file name in which %J is replaced by the job ID
#BSUB -e ACCBaseCode1.err  # error file name in which %J is replaced by the job ID
#BSUB -q gpgpu               # queue

cd /glade/p/work/pkommera/finalACCOptBaseCode1/homme-depot/src/full
#run the executable
mpirun.lsf ./dg2d
