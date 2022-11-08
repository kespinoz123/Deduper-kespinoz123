#!/bin/bash

#SBATCH --account=bgmp                  ### SLURM account which will be charged for the job
#SBATCH --job-name=Deduper              ### Job Name
#SBATCH --output=Deduper_%j.out         ### File in which to store job output
#SBATCH --error=Deduper-%j.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00               ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                       ### Node count required for the job
#SBATCH --cpus-per-task=4               ### Number of cpus (cores) per task
#SBATCH --partition=bgmp                ### partition to run things

in_sam="/projects/bgmp/kespinoz/bioinfo/Bi624/Deduper/C1_SE_uniqAlign.sam"
umi="/projects/bgmp/kespinoz/bioinfo/Bi624/Deduper/STL96.txt"
out_sam="/projects/bgmp/kespinoz/bioinfo/Bi624/Deduper/C1_SE_uniqAlign.deduped.sam"

/usr/bin/time -v ./Espinoza_deduper.py -f $in_sam -u $umi -o $out_sam