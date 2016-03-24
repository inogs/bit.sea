#! /bin/bash


LOC=/gpfs/work/Isc\rC_MYMEDBIO/COPERNICUS/online_validation_data
INPUTDIR=$LOC/output_phys
 ORIGDIR=$LOC/LINKED_ORIG/
GEN__DIR=$LOC/FORCING_GEN
 OUT_DIR=$LOC/FORCINGS/

mkdir -p $ORIGDIR
mkdir -p $GEN__DIR
mkdir -p $OUT_DIR

sed -e "s%@@ORIGDIR@@%${ORIGDIR}%g"  -e "s%@@OUTDIR@@%${OUT_DIR}%g" namelist.tpl > $GEN__DIR/namelist

cd $ORIGDIR
for I in `ls $INPUTDIR/*nc`; do
   filename=`basename $I` 
   new_name=${filename:21:6}_${filename:30:4}
   ln -fs $I $new_name
done




cd $ORIGDIR

ls *U.nc > $GEN__DIR/nomefile_U
ls *V.nc > $GEN__DIR/nomefile_V
ls *W.nc > $GEN__DIR/nomefile_W
ls *T.nc > $GEN__DIR/nomefile_T

cd $GEN__DIR
cp /pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/HOST/pico/bin/ForcingGenerator .

module purge
module load autoload intelmpi/5.0.1--binary mkl/11.2.0--binary
ulimit -s unlimited
./ForcingGenerator

