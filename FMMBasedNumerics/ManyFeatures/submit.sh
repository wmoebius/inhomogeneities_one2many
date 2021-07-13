#!/bin/bash

# start jobs
for rhoid in $(seq 0 $(python -c "import param; print(param.rhos.size-1)"))
do
    for gammaid in $(seq 0 $(python -c "import param; print(param.gammas.size-1)"))
    do
        cat runonodyssey.slurm | sed s/RHOID/${rhoid}/g | sed s/GAMMAID/${gammaid}/g > runnow.slurm
        sbatch --no-requeue --array=1-64 runnow.slurm
        rm runnow.slurm
        sleep 7
    done
done
