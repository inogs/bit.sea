GSSOBS_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/
annaDIR=/pico/scratch/userexternal/ateruzzi/Sat_Check_Clima/

DIR_ORIG1km=$GSSOBS_DIR/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/ORIG/
DIR_CLIMA1km=$GSSOBS_DIR/CLIMATOLOGY/SAT/CCI_1km/
FILEClima=$DIR_CLIMA1km/SatClimatology.nc

MESH=SAT1km_mesh

## Creation of subask indexes for sat mesh
# To be executed only one time (once for each climatology)

# OUT_SUBMASK=$DIR_CLIMA1km 
##### Io metterei là questi file perché sonno basati sulla climatologia

OUT_SUBMASK=$annaDIR/SUBMASKsat
# mkdir -p $OUT_SUBMASK

# echo python sat_indsub.py -o $OUT_SUBMASK -c $FILEClima -m $MESH

# --------------------------------
## Check of sat files

OUT_CHECK=$annaDIR/OUTPUT # 
OUT_STATS=$OUT_CHECK/STATISTICS
mkdir -p $OUT_CHECK/CHECKED
mkdir -p $OUT_CHECK/REJECTED
mkdir -p $OUT_STATS

echo python sat_check.py -i $DIR_ORIG1km -o $OUT_CHECK -c $FILEClima -m $MESH  -s $OUT_SUBMASK -w $OUT_STATS

# --------------------------------
## Figures

# OUT_FIG=FIGURES/$TYPEIN/
# mkdir -p $OUT_FIG
# echo python plot_climastats_time.py -o $OUT_FIG -i $OUT_STATS -m $OUT_SUBMASK


# --------------------------------
## Ave sat

CHECKMULTI=$GSSOBS_DIR/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/CHECKED/
WEEKLYMULTI=$annaDIR/WEEKLY
mkdir -p $WEEKLYMULTI/AVEfiles/
mkdir -p $WEEKLYMULTI/AVEdates/
echo python aveSat.py -i $CHECKMULTI -o $WEEKLYMULTI -m SAT1km_mesh -t weekly_tuesday


# --------------------------------
## Compose blacklisting file

DADIR=$annaDIR/FakeArchive
OUTNC=$annaDIR/NC_OUT
MASKMOD=/pico/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
mkdir -p $OUTNC

echo python blacklisting_ncfile.py -r $OUT_CHECK/REJECTED -d $DADIR -s $WEEKLYMULTI/AVEdates -o $OUTNC -t SAT1km_mesh -m $MASKMOD

