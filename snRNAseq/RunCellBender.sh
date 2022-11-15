#! /bin/bash

#$ -cwd
#$ -e ErrFiles/CB.$TASK_ID.err
#$ -o ErrFiles/CB.$TASK_ID.log
#$ -l h_vmem=20g
#$ -l h_rt=120:00:00
#$ -l os="RedHat7"
#$ -pe smp 8
#$ -binding linear:8
#$ -R y
#$ -t 1-numReps

input=samp_$SGE_TASK_ID
outdir=$input/CellBender
mkdir $outdir
numCells=10000

cellbender remove-background --input $input --output $outdir/cellbender_out.h5 --expected-cells $numCells --total-droplets-included 25000