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

#baseppcon--> da dove lancio 
BASEPPCON=/g100_scratch/userexternal/camadio0/PPCON/CODE_loss_attention_max_PPCon_CA/
INDIR=$BASEPPCON/ds/
# Creo il file tensor csv
python -u clustering/clustering.py -i $INDIR 
# clustering.py --> utils_clustering.py --> clustering/make_ds_clustering.py


# appendo il dataset
#qui o sovrascrivo SUPERFLOAT o
#creo una copia su cui lavoro e nomino $SUPERFLOAT_DIR
INDIR=$BASEPPCON/ds/ # come sopra ma la scrivo per chiarezza
metadata_OUTDIR=$BASEPPCON/results/
SUPERFLOAT_DIR=SUPERFLOAT -->  overwrite  the input dataset
mkdir -p $SUPERFLOAT_DIR 
python -u make_generated_ds/generate_netcdf_netcdf4.py -i $INDIR -o $metadata_OUTDIR -s $SUPERFLOAT_DIR


#
#python superfloat_chla.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
#python superfloat_o2o.py -Tst $DATE_start -Tend $DATE_end -Tday $DATE_DAY -o $OUTDIR -v $VARNAME
#python superfloat_nitrate.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
#python superfloat_par.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
#python superfloat_bbp700.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
#python superfloat_ph.py -s $DATE_DAY -e $DATE_end -o $OUTDIR -f
#python dump_index.py -i $OUTDIR -o ${OUTDIR}/Float_Index.txt -t superfloat
#
