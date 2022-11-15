#! /bin/bash -l

#$ -cwd
#$ -l h_vmem=50g
#$ -e aggr3.0.err
#$ -o aggr3.0.out
#$ -l h_rt=18:00:00
#$ -l os=RedHat7

cellranger aggr --id=EXPR --csv=/path/to/aggr_sampleSheet.csv --normalize=mapped