. ../online/profile.inc

export MASKFILE=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc
#export MASKFILE="/g100_scratch/userexternal/camadio0/" + run + "/wrkdir/MODEL/meshmask.nc" # It is the same of that one in the previous line

run=MedBFM4.2_Q24_v3 


TABLES_DIR=/g100_work/OGS_devC/Benchmark/pub/lfeudale/HPLC/Test_Benchmark/$run/tables/
FIGURES_DIR=/g100_work/OGS_devC/Benchmark/pub/lfeudale/HPLC/Test_Benchmark/$run/figures/
STAT_PROFILES_DIR=/g100_scratch/userexternal/lfeudale/dev/TRANSITION_V11/STATIC/STAT_PROFILES_MedBFM4.2_Q24_v3/

# CREATE Directories for PFTS validation with HPLC:
opa_prex_or_die "mkdir -p $TABLES_DIR $FIGURES_DIR "

#
# GENERATES THE ".csv" FILES FOR MODEL AND REF DATA WITH MEAN AND STD VALUES FOR MED BASIN
# FOR SUMMER AND WINTER:
opa_prex_or_die " python static_clim_validation_HPLC.py -i $STAT_PROFILES_DIR -o $TABLES_DIR -m $MASKFILE -s 20190101 -e 20200101"

# plot the results for every pfts:

opa_prex_or_die "python plot_HPLC_clim.py -i $TABLES_DIR -o $FIGURES_DIR -p P1l"
opa_prex_or_die "python plot_HPLC_clim.py -i $TABLES_DIR -o $FIGURES_DIR -p P2l"
opa_prex_or_die "python plot_HPLC_clim.py -i $TABLES_DIR -o $FIGURES_DIR -p P3l"
opa_prex_or_die "python plot_HPLC_clim.py -i $TABLES_DIR -o $FIGURES_DIR -p P4l"
