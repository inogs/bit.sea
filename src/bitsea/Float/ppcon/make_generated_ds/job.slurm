#! /bin/bash

#SBATCH --job-name=generate_superfloat
#SBATCH --ntasks=1
#SBATCH --time=01:05:00
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_interactive
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=camadio@inogs.it

date
python -u generate_netcdf_netcdf4.py

date
