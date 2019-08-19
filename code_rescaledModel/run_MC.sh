#!/bin/sh
#
#SBATCH --job-name=th150N10
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=60GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dq271@nyu.edu


cd /home/dq271/tKdV_simu/
matlab -nodisplay -r "SolverKdV_SymplecticM4a_MC(1E4, 0.,10, -150.,1., 1)"
