### Float

RUN=V7C/TEST_07/

SCRATCHDIR=/gpfs/scratch/userexternal/gbolzon0/

      DA_DIR=$SCRATCHDIR/$RUN/wrkdir/MODEL/DA__FREQ_1/
  PREPROC_DIR=$SCRATCHDIR/$RUN/wrkdir/DA/
  #PREPROC_DIR=/gpfs/work/OGS_prod_0/OPA/V6C/devel/wrkdir/analysis/2/DA
POSTFLOAT_DIR=/gpfs/scratch/userexternal/ateruzzi/$RUN/wrkdir/postproc_float/
      OUT_DIR=$POSTFLOAT_DIR/OUTSTATS/
     #OUT_RMSD=$POSTFLOAT_DIR/OUTRMSD/

mkdir -p $OUT_DIR

# cd $POSTFLOAT_DIR


echo python post_check_float.py -i $PREPROC_DIR -d $DA_DIR -o $OUT_DIR


OUT_FIG=$POSTFLOAT_DIR/OUTFIG/
mkdir -p $OUT_FIG

echo python plot_post_check_float.py -i $OUT_DIR -o $OUT_FIG


### Sat

SATSTATS_DIR=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V6C/SAT/MULTISENSOR/1Km/NRT/DAILY/STATISTICS
SUBMASKDIR=/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI_1km/SUBMASKsat/
OUTDIR=/gpfs/scratch/userexternal/ateruzzi/$RUN/wrkdir/postproc_sat/OUTSTATS
OUTFIG=/gpfs/scratch/userexternal/ateruzzi/$RUN/wrkdir/postproc_sat/OUTFIG
MASKFILE=$SCRATCHDIR/$RUN/wrkdir/MODEL/meshmask.nc

mkdir -p $OUTDIR
mkdir -p $OUTFIG

echo python post_check_sat.py -d $DA_DIR -o $OUTDIR -m $MASKFILE

echo python plot_check_sat.py -i $OUTDIR -o $OUTFIG -s $SATSTATS_DIR -b $SUBMASKDIR

