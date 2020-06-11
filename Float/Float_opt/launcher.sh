#! /bin/bash

#SBATCH --job-name=POST
#SBATCH -N1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=3gb
#SBATCH --account=OGS_dev_1
#SBATCH --partition=gll_usr_prod

cd $SLURM_SUBMIT_DIR

module purge
module load profile/base
module load intel/pe-xe-2018--binary intelmpi/2018--binary
module load autoload
module load hdf5/1.8.18--intel--pe-xe-2018--binary netcdf/4.6.1--intel--pe-xe-2018--binary
module load mpi4py/3.0.0--intelmpi--2018--binary
source /gpfs/work/OGS20_PRACE_P/COPERNICUS/py_env_2.7.12/bin/activate
export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

INPUTDIR=/gpfs/scratch/userexternal/eterzic0/BGC-ARGO-DATA
OUTDIR=/gpfs/scratch/userexternal/gbolzon0/plazzari/Float_opt_2020
FLOAT_INDEX=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/FLOAT_BIO/Float_Index.txt

echo python Float_opt_converter_2.py -i $INPUTDIR -o $OUTDIR -f $FLOAT_INDEX


exit 0
# caso Float_opt_2019
INPUTDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT_2019/history/
OUTDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT_2019/
FLOAT_INDEX=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/FLOAT_LOVBIO/Float_Index.txt

python Float_opt_converter_2.py -i $INPUTDIR -o $OUTDIR -f $FLOAT_INDEX