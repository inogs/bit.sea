GSSOBS_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/

DIR_ORIG1km=$GSSOBS_DIR/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/ORIG/
DIR_CLIMA1km=$GSSOBS_DIR/CLIMATOLOGY/SAT/CCI_1km/
FILEClima=$DIR_CLIMA1km/SatClimatology.nc

MESH=SAT1km_mesh

OUT_SUBMASK=SubMASK/
mkdir -p $OUT_SUBMASK

echo python sat_indsub.py -o $OUT_SUBMASK -c $FILEClima -m $MESH

echo -----------

OUT_DIR=OUTPUT
OUT_STATS=STATS
mkdir -p $OUT_DIR/CHECKED
mkdir -p $OUT_DIR/REJECTED
mkdir -p $OUT_STATS

echo python sat_climastats_time.py -i $DIR_ORIG1km -c $FILEClima -m $MESH  -s $OUT_SUBMASK -o $OUT_DIR -w $OUT_STATS

echo -----------

OUT_FIG=FIGURES/$TYPEIN/
mkdir -p $OUT_FIG
echo python plot_climastats_time.py -o $OUT_FIG -i $OUT_STATS -m $OUT_SUBMASK

echo -----------

CHECKMULTI=/gss/gss_work/DRES_OGS_BiGe/Observations//TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/CHECKED/
WEEKLYMULTI=WEEKLY
mkdir $WEEKLYMULTI
echo python aveSat16_dates.py -i $CHECKMULTI -o $WEEKLYMULTI


echo ------------

DADIR=FakeArchive
OUTNC=NC_OUT
MASKMOD=/pico/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
mkdir -p $OUTNC

echo python blacklisting_ncfile.py -r $OUT_DIR/REJECTED -d $DADIR -s $WEEKLYMULTI -o $OUTNC -t SAT1km_mesh -m $MASKMOD

