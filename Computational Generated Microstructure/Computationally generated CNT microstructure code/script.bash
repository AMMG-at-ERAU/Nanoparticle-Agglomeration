#PBS -S /bin/bash

#PBS -l walltime=23:00:00
#PBS -l procs=37
#PBS -N matlab_test_CNT_Cond_parallel_1
#PBS -q normalq
#PBS -m abe
#PBS -o outCond2.txt
#PBS -e errCond2.txt
#PBS -M motagis@my.erau.edu 
#PBS -d /scratch/motagis/code

module load matlab/2016b


cd $PBS_O_WORKDIR
module load matlab/2016b
matlab  -r Samarth_Cond_CNT_2D_Vega_one_Vf  &> outCond13..log




