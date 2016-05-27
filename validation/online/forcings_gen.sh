#! /bin/bash


if [ $# -lt 2 ] ; then
  usage
  exit 1
fi

for I in 1; do
   case $1 in
      "-d" ) LOC=$2;;
        *  ) echo "Unrecognized option $1." ; usage;  exit 1;;
   esac
   shift 2
done



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
   new_name=${filename:20:6}_${filename:29:4}
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
module load profile/base
module load autoload intelmpi/5.0.1--binary mkl/11.2.0--binary
export NETCDF_INC=/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/HOST/pico/include
export NETCDF_LIB=/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/HOST/pico/lib
export LD_LIBRARY_PATH=/pico/home/usera07ogs/a07ogs00/OPA/V2C-dev/HOST/pico/lib:$LD_LIBRARY_PATH
ulimit -s unlimited
./ForcingGenerator
