#!/bin/bash
#SBATCH -J runwrapper #A single job name for the array
#SBATCH -o LOGFILES/runwrapper_rhoid_RHOID_gammaid_GAMMAID_job_%A_%a.out #Standard output
#SBATCH -e LOGFILES/runwrapper_rhoid_RHOID_gammaid_GAMMAID_job_%A_%a.err #Standard error
#SBATCH -p shared #Partition
#SBATCH -t 7-00:00 #Running time
#SBATCH --mem 5000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --open-mode=append
#SBATCH --mail-type=END,FAIL,REQUEUE      # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=w.moebius@exeter.ac.uk  #Email to which notifications will be sent

# get random number seed from file 
RANDNUMSEED=$(head -n ${SLURM_ARRAY_TASK_ID} fixedrandnumseeds_rhoid_RHOID_gammid_GAMMAID.txt | tail -n 1 | awk '{print $2}')

# https://www.shell-tips.com/2010/06/14/performing-math-calculation-in-bash/
# make sure not all tasks start at the same time
# vague recollection: jobs starting at the same time caused problems
sleep $((SLURM_ARRAY_TASK_ID*71))

# logging version of CentOS
cat /etc/centos-release

module purge
module load ffmpeg/4.0.2-fasrc01
module load Anaconda3/2020.11
__conda_setup="$('/n/sw/eb/apps/centos7/Anaconda3/2020.11/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate skfmm

echo "FMM started."
python fmmmanyfeatures.py RHOID GAMMAID $SLURM_ARRAY_TASK_ID $RANDNUMSEED
echo "FMM ended."
if [ $SLURM_ARRAY_TASK_ID -le 3 ]
then
  echo "Making video started."
  python make_video.py RHOID GAMMAID $SLURM_ARRAY_TASK_ID
  echo "Making video ended."
fi

# create command used for stating amount of memory used
# https://rc.fas.harvard.edu/resources/odyssey-quickstart-guide/
# https://slurm.schedmd.com/job_array.html
if [ $SLURM_ARRAY_TASK_ID -eq $SLURM_ARRAY_TASK_MIN ]
then
  echo "sacct -j $SLURM_ARRAY_JOB_ID --format=JobID,JobName,ReqMem,MaxRSS,Elapsed > summary_${SLURM_ARRAY_JOB_ID}.log" >> LOGFILES/makesummary.sh
fi
