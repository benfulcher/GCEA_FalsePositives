#!/bin/csh
#PBS -N mouseBrainNullComputeSBPspatial
#PBS -o mouseBrainNull_SBPspatial.txt
#PBS -q yossarian
#PBS -l nodes=1:ppn=12
#PBS -l mem=144GB
# Minimum acceptable walltime: day-hours:minutes:seconds
#PBS -l walltime=140:00:00
# Email user if job ends or aborts
#PBS -m ea
#PBS -M ben.fulcher@sydney.edu.au
#PBS -j oe
#PBS -V

# Show the host on which the job ran
hostname
module load Matlab2018a

# Move
cd $PBS_O_WORKDIR
cd ../../

# Launch the Matlab job:
set jobText = "startup; parpool('local',10);params = GiveMeDefaultParams('mouse','all');params.e.whatEnsemble = 'customEnsemble';NullComputation(params); exit"
matlab -nodesktop -r "$jobText"
