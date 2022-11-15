#! /bin/bash -l

#$ -cwd
#$ -e /path/to/output/directory/mkfastq/mkfastq3.0.err
#$ -o /path/to/output/directory/mkfastq/mkfastq3.0.out
#$ -l h_vmem=50g
#$ -l h_rt=18:00:00
#$ -l os=RedHat7

cellranger mkfastq --run=/path/to/rawData --samplesheet=/path/to/sampleSheet.csv 
 --output-dir=/path/to/output/directory/mkfastq/ --localmem=96 --jobmode=local
