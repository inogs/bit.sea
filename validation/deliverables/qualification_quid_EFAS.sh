# MATCHUP VALIDATION using EFAS data: taken from Descriptor of CMEMS-Med-biogeochemistry-ScQP-...REAN....docx

# QUID REANALYSIS
# SECTION 4:
export MASKFILE=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/meshmask.nc
export ONLINE_REPO=gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V10C/


######################### MATCHUP PUNTUALE ####################

OUTDIR=/g100_work/OGS_devC/Benchmark/pub/lfeudale/MedBFM4.2_Q24_v3/STATIC/COASTAL/matchup_validation/density_plots
#OUTDIR=matchup_validation_run04_SEASONS/density_plots
mkdir -p $OUTDIR/Dens_N1p $OUTDIR/Dens_N3n $OUTDIR/Dens_N4n $OUTDIR/Dens_N5s $OUTDIR/Dens_O2o $OUTDIR/Dens_P_l
ln -sf profiler_RA_N.py profiler_RA.py
echo python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_N1p  -v N1p -m 0 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_N1p  -v N1p -m 0 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_N3n  -v N3n -m 0 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_N4n  -v N4n -m 0 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_N5s  -v N5s -m 0 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_O2o  -v O2o -m 140 -c supercoastal
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_P_l  -v P_l -m 0 -c supercoastal

if [ 1 == 0 ]; then
mkdir -p $OUTDIR/Dens_DIC $OUTDIR/Dens_ALK $OUTDIR/Dens_pH $OUTDIR/Dens_pCO2
ln -sf profiler_RA_C.py profiler_RA.py
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_DIC -v DIC -m 2150 #1500 #2050
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_ALK -v ALK -m 2350
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_pH -v pH -m 7.9
python density_plots.py   -M $MASKFILE  -o $OUTDIR/Dens_pCO2 -v pCO2 -m 300
fi

OUTDIR=/g100_work/OGS_devC/Benchmark/pub/lfeudale/MedBFM4.2_Q24_v3/STATIC/COASTAL/matchup_validation/
#mkdir -p $OUTDIR/validation_tables/OpenSea $OUTDIR/validation_tables/Coast $OUTDIR/validation_tables/SuperCoastal/
mkdir -p $OUTDIR/validation_tables/SuperCoastal/
#python matchup_validation.py -o $OUTDIR/validation_tables/OpenSea  -m $MASKFILE -a "OpenSea"
#python matchup_validation.py -o $OUTDIR/validation_tables/Coast  -m $MASKFILE -a "Coast"
python matchup_validation.py -o $OUTDIR/validation_tables/SuperCoastal/ -m $MASKFILE -a "SuperCoastal"

#python matchup_validation.py -o /g100_scratch/userexternal/lfeudale/NUTRIENT_CARBON_newDATA/bit.sea2/bit.sea/validation/deliverables/matchup_validation_SUPERCOASTAL/density_plots/validation_tables/SuperCoastal/ -m $MASKFILE -a "SuperCoastal"


# CREATE PLOTS OF VERTICAL PROFILE MATCHUPS
OUTDIR=/g100_work/OGS_devC/Benchmark/pub/lfeudale/MedBFM4.2_Q24_v3/STATIC/COASTAL/matchup_validation/ #matchup_validation_run04_SEASONS
mkdir -p $OUTDIR/Vert_Prof/N1p $OUTDIR/Vert_Prof/N3p $OUTDIR/Vert_Prof/N4n $OUTDIR/Vert_Prof/N5s $OUTDIR/Vert_Prof/O2o $OUTDIR/Vert_Prof/P_l $OUTDIR/Vert_Prof/DIC $OUTDIR/Vert_Prof/ALK $OUTDIR/Vert_Prof/pH $OUTDIR/Vert_Prof/pCO2 
ln -sf profiler_RA_N.py profiler_RA.py
# vertical_profiles.py is only for OPEN_SEA AREAS!!! It creates plot of vertical profiles
echo python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/N1p  -v N1p -c supercoastal | grep corr > $OUTDIR/table.N1p_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/N1p  -v N1p -c supercoastal | grep corr > $OUTDIR/table.N1p_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/N3n  -v N3n -c supercoastal | grep corr > $OUTDIR/table.N3n_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/N4n  -v N4n -c supercoastal | grep corr > $OUTDIR/table.N4n_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/N5s  -v N5s -c supercoastal | grep corr > $OUTDIR/table.N5s_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/O2o  -v O2o -c supercoastal | grep corr > $OUTDIR/table.O2o_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/P_l  -v P_l -c supercoastal | grep corr > $OUTDIR/table.P_l_corr.txt

if [ 1 == 0 ]; then
ln -sf profiler_RA_C.py profiler_RA.py
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/DIC  -v DIC | grep corr > $OUTDIR/table.DIC_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/ALK  -v ALK | grep corr > $OUTDIR/table.ALK_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/pH   -v pH  | grep corr > $OUTDIR/table.pH_corr.txt
python vertical_profiles.py -m $MASKFILE -o $OUTDIR/Vert_Prof/pCO2 -v pCO2 | grep corr > $OUTDIR/table.pCO2_corr.txt
fi

exit 0
