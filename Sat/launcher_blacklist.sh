#! /bin/bash
annaDIR=/pico/scratch/userexternal/ateruzzi/Sat_Check_Clima/

  DIR_ORIG1km=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/ORIG/
DIR_CHECK_1km=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MULTISENSOR/1Km/NRT/DAILY/
DIR_CLIMA1km=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/
   FILEClima=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SatClimatology.nc
 DIR_SUBMASK=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SUBMASKsat

#DIR_CHECK_1km=/marconi_scratch/userexternal/gbolzon0/SAT_CHECK

## Creation of subask indexes for sat mesh
# To be executed once and for all (once for each climatology)
# echo python sat_indsub.py -o $DIR_SUBMASK -c $FILEClima -m SAT1km_mesh

# --------------------------------
## Creation of climatology statistics for the Mediterranean subbasins
# To be executed once and for all (once for each climatology)
# OUT_STATSCLIM=$annaDIR/STATS_CLIM/
# mkdir $OUT_STATSCLIM
# echo python stats_sub_clima.py -c $FILEClima -m SAT1km_mesh  -s $DIR_SUBMASK -w $OUT_STATSCLIM

# --------------------------------
## Check of sat files

OUT_STATS=$DIR_CHECK_1km/STATISTICS
mkdir -p $DIR_CHECK_1km/CHECKED
mkdir -p $DIR_CHECK_1km/REJECTED
mkdir -p $OUT_STATS

echo python sat_check.py -i $DIR_ORIG1km -o $DIR_CHECK_1km -c $FILEClima -m SAT1km_mesh  -s $DIR_SUBMASK -w $OUT_STATS

# --------------------------------
## Figures

# OUT_FIG=FIGURES/$TYPEIN/
# mkdir -p $OUT_FIG
# echo python plot_climastats_time.py -o $OUT_FIG -i $OUT_STATS -m $DIR_SUBMASK


# --------------------------------
## Ave sat

WEEKLYMULTI=$annaDIR/WEEKLY
WEEKLYDATES=$annaDIR/WEEKLY_2_AVEDATES
mkdir -p $WEEKLYMULTI
mkdir -p $WEEKLYDATES
echo python aveSat.py -i $DIR_CHECK_1km/CHECKED -o $WEEKLYMULTI -d $WEEKLYDATES -m SAT1km_mesh -t weekly_tuesday


# --------------------------------
## Compose blacklisting file

FLAGSATDIR=$annaDIR/FakeArchive
OUTNC=$annaDIR/NC_OUT
MASKMOD=/pico/scratch/userexternal/ateruzzi/MASKS16/meshmask.nc
mkdir -p $OUTNC

echo python blacklisting_ncfile.py -r $DIR_CHECK_1km/REJECTED -d $FLAGSATDIR -s $WEEKLYDATES -o $OUTNC -t SAT1km_mesh -m $MASKMOD

