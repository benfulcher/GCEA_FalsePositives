#!/bin/bash
# Set name of job shown in squeue
#SBATCH --job-name CFPR_humanSpatialLag
# Set project code account
#SBATCH --account=rn29
# Request CPU resources
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
# Memory usage (MB)
#SBATCH --mem-per-cpu=12000
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=120:00:00
# Email user if job fails or ends
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=ben.fulcher@sydney.edu.au
# Specify a queue (called a partition on SLURM)
# SBATCH --partition=m3a

# Set environment variables to run Matlab
module purge
module load matlab/r2018a

# Show the host on which the job ran and return to home repository directory
hostname
cd ../

# Show what SLURM ennvironment variables our environment has
env | grep SLURM

# Launch the Matlab job
matlab -nodesktop -r "startup; parpool('local',16);
                        params = GiveMeDefaultParams('human');
                        params.g.whatSurrogate = 'spatialLag';
                        params.nulls.customShuffle = 'none';
                        SurrogateEnrichment(params); exit"
