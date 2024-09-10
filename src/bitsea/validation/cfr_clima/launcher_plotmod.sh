gssDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/

DIRSTATSClima=$gssDIR/CCI_1km/SUBBASIN_STATISTICS/
# Directory containing climatology statistics for satellite observations
# (the same climatology used for check on satellite observations before DA)

#DIRSatDA=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/SAT_CORRECTED/
DIRSatDA=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/CFR_CLIMA/SAT_chl/
# Directory of assimilated satellite

# DIRCHLMOD=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/AVE_FREQ_1/
#DIRCHLMOD=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/DA__FREQ_1/
DIRCHLMOD=/gpfs/scratch/userexternal/ateruzzi/CHECK_DA/CFR_CLIMA/DA_chllinkall/
# Directory of chl output

MASKFILE=/gpfs/scratch/userexternal/ateruzzi/RA_COAST2017/wrkdir/MODEL/meshmask.nc

OUT_FIG=/gpfs/scratch/userexternal/ateruzzi/ELAB_RACOAST2017/CFR_CLIMAwith2016/
mkdir -p $OUT_FIG
echo python plot_climamod_time.py -o $OUT_FIG -c $DIRSTATSClima -s $DIRSatDA -i $DIRCHLMOD -m $MASKFILE #-t 3
# -t is optional and indicates the number of months used for the plot
# if not provided, all the available dates are used

echo -----------

