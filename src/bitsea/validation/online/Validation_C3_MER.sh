# SST:

cd $OPA_BITSEA/validation/online
# SAT DIR for Temperature:
SAT_DAILY_DIR=${ONLINE_REPO}/SAT/SST/NRT/DAILY/CHECKED_24
# Model dir: ave_files
inputdir=${ONLINE_VALIDATION_DIR}/ARCHIVE_FC_FOR_SAT/
# Output dir: name of the day to validate
outputdir=${ONLINE_VALIDATION_DIR}/SST

# OPA_RUNDATE is the actual date in which I am running the model
opa_prex_or_die "python RunSatValidation_MER.py -i $inputdir  -o $outputdir -s $SAT_DAILY_DIR  -c null  -m $MASKFILE  -r $OPA_RUNDATE -v votemper"
INPDIR_SAT_V_DAILY=${OPA_INPDIR_ROOT}/analysis/VALIDATION/SAT/SST/NRT/
opa_cp "$outputdir/Validation_*nc $INPDIR_SAT_V_DAILY"

##################
#CHL:

SAT_DAILY_DIR=${ONLINE_REPO}/SAT/CHL/NRT/DAILY/CHECKED_24

inputdir=${ONLINE_VALIDATION_DIR}/ARCHIVE_FC_FOR_SAT/
outputdir=${ONLINE_VALIDATION_DIR}/CHL
# VERIFICARE il CLIMDIR: temporary we set null
#climdir=${ONLINE_VALIDATION_DIR}/CHL_CLIM/2D/

#opa_prex_or_die "python RunSatValidation.py -i $inputdir  -o $outputdir -s $SAT_DAILY_DIR  -c $climdir  -m $MASKFILE  -r $OPA_RUNDATE -v P_l"
opa_prex_or_die "python RunSatValidation_MER.py -i $inputdir  -o $outputdir -s $SAT_DAILY_DIR  -c null  -m MASKFILE  -r $OPA_RUNDATE -v P_l"


INPDIR_SAT_V_DAILY=${OPA_INPDIR_ROOT}/analysis/VALIDATION/SAT/CHL/NRT/
opa_cp "$outputdir/Validation_*nc $INPDIR_SAT_V_DAILY"
opa_prex "python sat_plotter.py  -i ${INPDIR_SAT_V_DAILY} -o ${outputdir}   -d $OPA_RUNDATE -v P_l"
opa_cp " ${outputdir}/*png $ONLINE_VALIDATION_DIR/MEDEAF_PUB/Sat_validation"


# MAP CREATION, for SST and CHLa:
# for the SATmap, use as basis the script:
# https://github.com/inogs/bit.sea/blob/master/src/bitsea/validation/deliverables/sat_model_RMSD_and_plot.py
# for Model, use as basis "build_layer_map.py"
#mit_prex_or_die "mpirun python $MEDEAF_DIR/build_layer_map.py -i $NETCDF_DIR -o $OUTPUTDIR -d $MIT_RUNDATE -m $MASKFILE -p $PLOTLISTFILE "


