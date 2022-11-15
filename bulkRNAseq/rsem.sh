#!/bin/bash -l

INPUTSAMPLES=(x y z)
SAMPLENAME="${INPUTSAMPLES[$SGE_TASK_ID - 1]}"

#$ -cwd
#$ -l h_vmem=100g
#$ -l h_rt=48:00:00
#$ -t 1-numReps

rsem-calculate-expression -p 8  --bowtie2 --paired-end --estimate-rspd --append-names --sort-bam-by-coordinate /path/to/data/$SAMPLENAME/$SAMPLENAME.unmapped.1.fastq.gz /path/to/data/$SAMPLENAME/$SAMPLENAME.unmapped.2.fastq.gz /path/to/reference /path/to/data/$SAMPLENAME/$SAMPLENAME
