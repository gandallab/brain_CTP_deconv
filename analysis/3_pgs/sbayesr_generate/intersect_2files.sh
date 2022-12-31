#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=2:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc
module load R/4.0.2

## substitute the command to run your code
## in the two lines below:
echo '/usr/bin/time -v hostname'
/usr/bin/time -v hostname

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

#==============================================================================

## If using arguments with flags
## while getopts f1:c1:h1:f2:c2:h2:m:t:o:n1:n2 option
## do
## case "${option}"
## 	in 
## 	f1) f1n=${OPTARG};;
## 	c1) c1=${OPTARG};;
## 	h1) h1=${OPTARG};;
## 	f2) file2=${OPTARG};;
## 	c2) c2=${OPTARG};;
## 	h2) h2=${OPTARG};;
## 	mc) mc=${OPTARG};;
## 	ot) out_type=${OPTARG};;
## 	od) outdir=${OPTARG};;
## 	n1) f1n=${OPTARG};;
## 	n2) file2n=${OPTARG};;
## esac
## done

## If using arguments without flags
## ${f1}=$1
## ${c1}=$2
## ${h1}=$3
## ${file2}=$4
## ${c2}=$5
## ${h2}=$6
## ${mc}=$7
## ${jt}=$8
## ${ot}=$9
## ${od}=$10
## ${n1}=$11
## ${n2}=$12

echo File 1 directory is ${f1}
echo Taking file1 column ${c1}
echo Interpreting file1 header: ${h1}
echo File2 directory ${f2}
echo Taking file2 column ${c2}
echo Interpreting file2 header: ${h2}
echo Joining files over column name ${mc}
echo Choosing to ${ot} output files
echo Output directory is ${od}
echo Output names are ${n1} ${n2}

# If using qsub, can go straight here
Rscript ${scripts}/intersect_2files.R \
${f1} \
${c1} \
${h1} \
${f2} \
${c2} \
${h2} \
${mc} \
${ot} \
${od} \
${n1} \
${n2}

#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
