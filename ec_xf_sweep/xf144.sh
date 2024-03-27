#!/bin/bash

#SBATCH --account=jcaps1
#SBATCH --job-name=xf144oc
#SBATCH --mail-user=jwake@umich.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=13-23:59:00
#SBATCH --partition=standard-oc
#SBATCH --output=%x-%j.log


CASEFOLDER=/home/jwake/ngaec/ec_xf_sweep_144_oc
JOBFOLDER=$ESCRATCH/$SLURM_JOB_NAME-$SLURM_JOB_ID


cp -R $CASEFOLDER $JOBFOLDER
cd $JOBFOLDER

touch jobinfo
echo "Job ID $SLURM_JOB_ID" >> jobinfo
echo "Job Tasks $SLURM_NTASKS" >> jobinfo
echo "Job Nodes $SLURM_NODELIST ($SLURM_NUM_NODES total)" >> jobinfo
echo "Job Task Distribution $SLURM_TASKS_PER_NODE" >> jobinfo
echo "Case Folder $CASEFOLDER" >> jobinfo
echo "Case Folder Contents:" >> jobinfo
for f in `ls -1 $CASEFOLDER`; do
    echo "\t$f" >> jobinfo
done

echo "[$(date +%Y-%m-%d\ %H:%M:%S)] Starting job..."

srun ./nga.dp.gnu.opt.mpi.exe -i input_144 &> out-$SLURM_JOB_ID

echo "[$(date +%Y-%m-%d\ %H:%M:%S)] Done"


