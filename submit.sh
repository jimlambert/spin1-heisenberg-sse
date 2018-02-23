#!/bin/dash
# Script for call run script to submit jobs to SLURM scheduler.
# First export variable for run table
export TABLE=runtable.dat

# Count number of jobs to submit
N_CASES=$(cat "$TABLE" | wc -l)

for i in i$(seq 0 1 $N_CASES)
do
  echo ${i}
  LINE=`sed -n ${i}p "$TABLE"`
  echo ${LINE}
  pwd
  eval "$LINE"
done
