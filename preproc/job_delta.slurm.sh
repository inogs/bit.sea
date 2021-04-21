#! /bin/bash

#SBATCH --job-name=test_txt
#SBATCH -N2
#SBATCH --ntasks-per-node=20
#SBATCH --time=02:00:00
#SBATCH --mem=100gb
#SBATCH --account=OGS20_PRACE_P
#SBATCH --partition=gll_usr_prod
#SBATCH --qos=gll_qos_dbg

cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary

source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

INPUTDIR=/gpfs/scratch/userexternal/gcoidess/FORCINGS_REA24/1993_transf/
OUTPUTDIR=/gpfs/scratch/userexternal/gcoidess/FORCINGS_REA24/TXT/
INGVMASK=/gpfs/work/IscrC_REBIOMED/REANALISI_24/PREPROC/MASK/ogstm/meshmask_INGVfor_ogstm.nc


mpirun -np 40 python delta_t_from_uvw.py -i $INPUTDIR -o $OUTPUTDIR -m $INGVMASK

