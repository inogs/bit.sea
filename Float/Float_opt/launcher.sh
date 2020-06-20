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
#OUTDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT_2020/
FLOAT_INDEX=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/FLOAT_BIO/Float_Index.txt

python Float_opt_converter_2.py -i $INPUTDIR -o $OUTDIR -f $FLOAT_INDEX
cd ..
# Float_Index.0.txt structure 6901772/SR6901772_139.nc,36.174896,18.768433,20200514-10:21:00, SALI TEMP BBP700
python dump_index.py -i $OUTDIR -o $OUTDIR/Float_Index.0.txt -t Float_opt_20

### optional -- positions are supposed to be good
# edit check_time_pos.py in order to have float_dataset="FLOAT_BIO
echo python check_time_pos.py -i $OUTDIR/Float_Index.0.txt -o $OUTDIR/Float_Index.txt



# superfloat generation
export STATIC_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C

cp -r $OUTDIR /gpfs/scratch/userexternal/gbolzon0/plazzari/SUPERFLOAT/
# T,S, BBP are "generated" by copying from Float_opt_2020

SUPERFLOAT_DIR=/gpfs/scratch/userexternal/gbolzon0/plazzari/SUPERFLOAT
COMMON="-s 20120101 -e 20200601 -o $SUPERFLOAT_DIR"
python superfloat_chla.py    $COMMON
python superfloat_irr380.py  $COMMON
python superfloat_irr412.py  $COMMON
python superfloat_irr490.py  $COMMON
python superfloat_par.py     $COMMON

python superfloat_oxygen.py  $COMMON
python superfloat_ph.py      $COMMON
python superfloat_cdom.py    $COMMON
python superfloat_nitrate.py $COMMON


python dump_index.py -i $SUPERFLOAT_DIR -o $SUPERFLOAT_DIR/Float_Index.0.txt -t static_superfloat
python check_time_pos.py -i $SUPERFLOAT_DIR/Float_Index.0.txt -o $SUPERFLOAT_DIR/Float_Index.txt # no changes

exit 0
# caso Float_opt_2019
INPUTDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT_2019/history/
OUTDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/Float_OPT_2019/
FLOAT_INDEX=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/FLOAT_LOVBIO/Float_Index.txt

python Float_opt_converter_2.py -i $INPUTDIR -o $OUTDIR -f $FLOAT_INDEX
