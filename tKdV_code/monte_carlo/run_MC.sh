#!/bin/sh
#
#SBATCH --job-name=J32g024MC100_L3D024E100u0
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=60GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dq271@nyu.edu


cd /home/dq271/tKdV_simu/simu_MC
matlab -nodisplay -r "SolverKdV_SymplecticM4a_MC(32, 3,.24,100, 0., 1E4,-.5,.24)"
