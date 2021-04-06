BASEDIR=/gpfs/scratch/userexternal/ateruzzi/REA_24_ELAB/VALIDATION/TEST14/deliverables/
INDIR=$BASEDIR/matchup_validation/density_plots/validation_tables/OpenSea/
OUTDIR=TEST14/mathcup_compact_tables/

mkdir -p $OUTDIR

VAR=N3n

echo python compact_matchup_tables.py -i $INDIR -o $OUTDIR -v $VAR

INDIR=$BASEDIR/matchup_validation/density_plots/validation_tables/Coast/
echo python compact_matchup_tables_coast.py -i $INDIR -o $OUTDIR -v $VAR