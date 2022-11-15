#! /bin/bash -l

#$ -cwd
#$ -e /path/to/output/directory/mkfastq/REPLICATE/cellranger-3.0.err
#$ -o /path/to/output/directory/mkfastq/REPLICATE/cellranger-3.0.out
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -pe smp 8
#$ -binding linear:8
#$ -R y
#$ -l os=RedHat7

cellranger count --id=REPLICATE --transcriptome=/path/to/reference/ --expect-cells=NUMCELLS --fastqs=/path/to/output/directory/mkfastq/REPLICATE --sample=REPLICATE --chemistry=SC3Pv3 --localmem=64
