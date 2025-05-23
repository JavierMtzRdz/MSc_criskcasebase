#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=256M
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="setting10"

for i in {1..200}

do
     ID=${i}
     RUNSCRIPTS=~/projects/def-gcohenfr/ntamvada/july5_beluga6to11/runscripts/setting10
     LOGS=~/projects/def-gcohenfr/ntamvada/july5_beluga6to11/logs/setting10
echo "#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=200GB
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=${LOGS}/${ID}.std_out.2.txt
#SBATCH --job-name="setting10"

module load StdEnv/2020 r/4.2.2

Rscript setting10.R

">$RUNSCRIPTS/${ID}.run.job
echo "Submitting $RUNSCRIPTS/${ID}.run.job to the cluster"
sbatch $RUNSCRIPTS/${ID}.run.job
done
