# Descriptor of CMEMS-Med-biogeochemistry-ScQP-v3.1.docx

# QUID REANALYSIS
# SECTION 4:
. ../online/profile.inc

export MASKFILE=/g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/MASK/meshmask.nc
export ONLINE_REPO=/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/

#         INPUTDIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/MODEL/AVE_FREQ_2/
         INPUTDIR=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/AVE_FREQ_2/
#   INPUT_AGGR_DIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP
   INPUT_AGGR_DIR=/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v3/wrkdir/MODEL/AVE_FREQ_2/
MASKFILE="/g100_scratch/userexternal/camadio0/MedBFM4.2_Q24_v2/wrkdir/MODEL/meshmask.nc"

STAT_PROFILES_DIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/

#SAT_KD_WEEKLY_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/KD490/DT/SIGMA_2.0/WEEKLY_4_24
SAT_KD_WEEKLY_DIR=/g100_scratch/userexternal/lfeudale/KD490/SAT/KD490/DT/SIGMA_2.0_KDmin0.022/WEEKLY_4_24/
SAT_CHLWEEKLY_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/CHL/DT/WEEKLY_4_24
SAT_CHLPFTsWEEKLY_DIR=/g100_scratch/userexternal/gbolzon0/V10C/SAT/CHL/DT/WEEKLY_4_24
KD_MODEL_DIR=/g100_scratch/userexternal/gbolzon0/V10C/run4.19/wrkdir/POSTPROC/output/AVE_FREQ_3/KD_WEEKLY/
BASEDIR=/g100_scratch/userexternal/lfeudale/validation/V10C/run4.19/PROFILATORE/


# CREATE Directories for SAT:
opa_prex_or_die "mkdir -p Fig4.2 Fig4.3/offshore Fig4.3/coast table4.1 table4.2"


for VAR in   kd490  P_l  P1l    P2l    P3l    P4l  ; do
   OPENSEA_PKL=${VAR}_open_sea.pkl
     COAST_PKL=${VAR}_coast.pkl

   MODELDIR=$INPUTDIR
   SAT_DIR=$SAT_CHLWEEKLY_DIR
   if [ $VAR == 'P_l' ] ; then 
        MODELDIR=$INPUT_AGGR_DIR #; fi
        SAT_DIR=$SAT_CHLWEEKLY_DIR
    fi
   if [ $VAR == 'kd490' ] ; then 
        MODELDIR=$KD_MODEL_DIR
        SAT_DIR=$SAT_KD_WEEKLY_DIR
    fi
 
   if [$VAR == P1l] | [$VAR == P2l] | [$VAR == P3l] | [$VAR == P4l]; then 
   # Plot PFTs comparison with HPLC dataset:
       INDIR_HPLC_CLIM=/g100_scratch/userexternal/lfeudale/validation/V10C/run4.19/bit.sea/validation/deliverables/
       opa_prex_or_die "python plot_HPLC_clim.py -i $INDIR_HPLC_CLIM -p $VAR "
    fi

   LAYER=10
   # if [ $VAR == 'kd490' ] ; then LAYER=0 ; fi
   
   
   opa_prex_or_die "python ScMYvalidation_plan.py -v $VAR -s $SAT_DIR -i $MODELDIR -m $MASKFILE -c open_sea -l $LAYER  -o $OPENSEA_PKL"
   opa_prex_or_die "python ScMYvalidation_plan.py -v $VAR -s $SAT_DIR -i $MODELDIR -m $MASKFILE -c coast    -l $LAYER  -o $COAST_PKL"
   opa_prex_or_die "python plot_timeseries_STD.py -v $VAR -i $OPENSEA_PKL -o ./Fig4.2/ "
   opa_prex_or_die "python plot_timeseries_RMS_CORR.py -v $VAR -i $OPENSEA_PKL -o Fig4.3/offshore " # table4.1
   opa_prex_or_die "python plot_timeseries_RMS_CORR.py -v $VAR -i $COAST_PKL -o Fig4.3/coast  "   # table4.2
   
   opa_prex_or_die "cp Fig4.3/offshore/table4.1_${VAR}.dat table4.1 "
   opa_prex_or_die "cp Fig4.3/coast/table4.1_${VAR}.dat    table4.2/table4.2_${VAR}.dat "

done

mkdir -p pfts/
opa_prex_or_die "python plot_pfts.py -i $PWD -o $PWD/pfts/ "


######################

# CREATE MAP PPN INTEGRAL 0-200m (python2 che python3 ha un errore sul name "long")
#opa_prex_or_die "mkdir -p Fig4.7refScale"
#opa_prex_or_die "python averager_and_plot_map_ppn_refScale.py -i $INPUT_AGGR_DIR  -v ppn  -t integral -m $MASKFILE -o Fig4.7refScale -l Plotlist_bio.xml -s 20190101 -e 20200101"

run="MedBFM4.2_Q24_v3"
OUTDIR="/g100_work/OGS_devC/Benchmark/pub/lfeudale/MedBFM4.2_Q24_v3/PPN"
# CREATE MAP PPN INTEGRAL 0-200m (python2 che python3 ha un errore sul name "long")
opa_prex_or_die "mkdir -p Fig4.7refScale"
#opa_prex_or_die "python averager_and_plot_map_ppn_refScale.py -i $INPUT_AGGR_DIR  -v ppn  -t integral -m $MASKFILE -o Fig4.7refScale -l Plotlist_bio.xml -s 20190101 -e 20200101"
opa_prex_or_die "python averager_and_plot_map_ppn_refScale_16basin_RMSD.py -i $INPUT_AGGR_DIR  -v ppn  -t integral -m $MASKFILE -o $OUTDIR -l Plotlist_bio.xml -s 20190101 -e 20200101"



# FLOAT SECTION:

# HOVMOELLER:
opa_prex_or_die "mkdir -p Fig4.4 Fig4.6 Fig4.12 Fig4.13 Fig4.15 tmp_nc table4.4 table4.8"
NCDIR=tmp_nc
OUTDIR=Fig4.4
opa_prex_or_die "echo python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR"
opa_prex_or_die "python SingleFloat_vs_Model_Stat_Timeseries.py -m $MASKFILE -b $BASEDIR -o $NCDIR"
opa_prex_or_die "echo python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR"
opa_prex_or_die "python Hov_Stat_plot.py -m $MASKFILE -i $NCDIR -o $OUTDIR -b $BASEDIR"
opa_prex_or_die "cp $OUTDIR/*N3n*.png Fig4.12"
opa_prex_or_die "cp $OUTDIR/*O2o*.png Fig4.15"

#------------------------------------------------------
# Statistics about key processes for different basins
#
# Figures 4.6 and 4.13 + tables 4.4 and 4.8
# CHL-PROF-D-CLASS4-PROF-CORR-BASIN
# NIT-PROF-D-CLASS4-PROF-CORR-BASIN
#  DO-PROF-D-CLASS4-PROF-CORR-BASIN 

OUTDIR=Fig4.6
#mkdir -p $NCDIR/chla_check/ $NCDIR/nitrate_check/
opa_prex_or_die "python BASIN_Float_vs_Model_Stat_Timeseries_monthly.py -m $MASKFILE -o $NCDIR"
opa_prex_or_die "python BASIN_Float_vs_Model_Stat_Timeseries_monthly_plotter.py -m $MASKFILE -i $NCDIR -o $OUTDIR"
opa_prex_or_die "cp $OUTDIR/N3n*.png Fig4.13"

opa_prex_or_die "cp $OUTDIR/P_l_tab_statistics_SHORT.txt table4.4/"
opa_prex_or_die "cp $OUTDIR/N3n_tab_statistics_SHORT.txt table4.8/"

#-----------------------------------------
# BIAS and RMSD averaged on BASIN:
#
# BIOFLOATS SECTION: statistics on layers
# CHL-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
# NIT-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN
#  DO-LAYER-D-CLASS4-PROF-[BIAS/RMS]-BASIN


OUTFIGDIR=Floats_bias_rmse_Timeseries     # 8layer x 7sub x 3var = 168 png files
TABLE_DIR=Floats_bias_rmse_tables         #: 2stats x 3var        = 6 txt files, TABLE.O2o_BIAS.txt  with time average for each layer,sub
opa_prex_or_die "mkdir -p $OUTFIGDIR $TABLE_DIR table4.3/ table4.9/ table4.12/ Fig4.5/ Fig4.14/ Fig4.16/"
opa_prex_or_die "python biofloats_ms.py  -m $MASKFILE -o float_bias_rmse.nc"
opa_prex_or_die "python biofloats_ms_plotter.py -i float_bias_rmse.nc -f $OUTFIGDIR -t $TABLE_DIR"
opa_prex_or_die "cp $TABLE_DIR/P_l_BIAS.txt $TABLE_DIR/P_l_RMSE.txt table4.3/"
opa_prex_or_die "cp $TABLE_DIR/N3n_BIAS.txt $TABLE_DIR/N3n_RMSE.txt table4.9/"
opa_prex_or_die "cp $TABLE_DIR/O2o_BIAS.txt $TABLE_DIR/O2o_RMSE.txt table4.12/"
opa_prex_or_die "cp $OUTFIGDIR/*P_l* Fig4.5"
opa_prex_or_die "cp $OUTFIGDIR/*N3n* Fig4.14"
opa_prex_or_die "cp $OUTFIGDIR/*O2o* Fig4.16"



####
# PLOTTING LAYER MAPS:

#OUTDIR=spaghettiplots
#mkdir -p $OUTDIR
#python plot_layer_timeseries_on_profiles.py -i $STAT_PROFILES_DIR -m $MASKFILE -o $OUTDIR

opa_prex_or_die "mkdir -p Fig4.1/ Fig4.7 Fig4.9 Fig4.10 Fig4.18 Fig4.19 Fig4.7bis"


COMMONS_PARAMS="-m $MASKFILE -o LayerMaps/  -l Plotlist_bio.xml -s 20190101 -e 20200101"

opa_prex_or_die "mkdir -p LayerMaps"
opa_prex_or_die "python averager_and_plot_map.py -i $INPUT_AGGR_DIR  -v P_l  -t mean $COMMONS_PARAMS"      # Fig4.1 CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v N3n  -t mean $COMMONS_PARAMS"     # Fig4.10 NIT-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v N1p  -t mean $COMMONS_PARAMS"      # Fig4.9  PHOS-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v N4n  -t mean $COMMONS_PARAMS"
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v N5s  -t mean $COMMONS_PARAMS"
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v O2o  -t mean $COMMONS_PARAMS"  # DO-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN --> Not showed for lack of ref
opa_prex_or_die "python averager_and_plot_map_ppn.py -i $INPUT_AGGR_DIR  -v P_c  -t mean -m $MASKFILE -o LayerMaps/  -l Plotlist_bio_Int.xml -s 20190101 -e 20200101"
opa_prex_or_die "python averager_and_plot_map_ppn.py -i $INPUT_AGGR_DIR  -v Z_c -t integral $COMMONS_PARAMS"

opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v ALK   -t mean $COMMONS_PARAMS"   # Ac-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN  --> not requested in the ScQP
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v DIC  -t mean $COMMONS_PARAMS"   # DIC-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN --> not requested in the ScQP
opa_prex_or_die "python averager_and_plot_map.py -i $INPUTDIR        -v CO2airflux -t mean $COMMONS_PARAMS"
opa_prex_or_die "cp LayerMaps/*P_l* Fig4.1/"
opa_prex_or_die "cp LayerMaps/*N3n* Fig4.10/"
opa_prex_or_die "cp LayerMaps/*N1p* Fig4.9/"
opa_prex_or_die "cp LayerMaps/*ppn* Fig4.7/"
opa_prex_or_die "cp LayerMaps/*ALK* Fig4.18/"
opa_prex_or_die "cp LayerMaps/*DIC* Fig4.19/"


#CHL-LAYER-Y-CLASS1-[CLIM/LIT]-MEAN from SATELLITE:
opa_prex_or_die "python sat_ave_and_plot.py      -i $SAT_CHLWEEKLY_DIR -m $MASKFILE  -o Fig4.1/"
opa_prex_or_die "python sat_model_RMSD_and_plot.py -s $SAT_CHLWEEKLY_DIR -i $INPUT_AGGR_DIR -m $MASKFILE  -o Fig4.1/ -st 20190101 -et 20200101"

#mkdir -p Fig4.8
#INTEGRALS_PPN=/g100_scratch/userexternal/gbolzon0/V9C/2019/TEST_01/wrkdir/POSTPROC/out/AVE_FREQ_2/INTEGRALS/PPN/
#python read_ppn_from_avescan_do_plot.py -c open_sea   -i $INTEGRALS_PPN -o Fig4.8


#########################   static dataset climatology section ###################################
# Fig. 4.11:
# PHO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
# NIT-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
#  DO-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers

# Fig. 4.19:
# ALK-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
# DIC-LAYER-Y-CLASS4-CLIM-BIAS/RMSD/CORR calculated on 8 layers
 
# Fig. 4.17:
# PH-PROF-Y-CLASS4-[CLIM/LIT]-MEAN
# pCO-PROF-Y-CLASS4-[CLIM/LIT]-MEAN

############
# COMPARISON WITH EMODnet CLIMATOLOGY:

#STAT_PROFILES_MONTHLY_DIR=/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v20_1/wrkdir/POSTPROC/output/MONTHLY/STAT_PROFILES/
#STAT_PROFILES_MONTHLY_DIR=/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/TEST_05/wrkdir/POSTPROC/output/AVE_FREQ_2/STAT_PROFILES/
# Figures 4.11 and 4.17 (previously named 4.19)
opa_prex_or_die "mkdir -p sim_vs_clim_profiles/ Fig4.11 Fig4.17 sim_vs_clim_profiles_OpenSea"
#python simulation_vs_clim.py -i $STAT_PROFILES_MONTHLY_DIR -o sim_vs_clim_profiles/ -s 20170101 -e 20180101 -m $MASKFILE
opa_prex_or_die "python simulation_vs_clim_extended.py  -i $STAT_PROFILES_DIR -o sim_vs_clim_profiles/ -s 20190101 -e 20200101 -m $MASKFILE"
opa_prex_or_die "mkdir -p sim_vs_clim_profiles_OpenSea"
opa_prex_or_die "python simulation_vs_clim_extended_OpenSea.py -i $STAT_PROFILES_DIR -o sim_vs_clim_profiles_OpenSea -s 20190101 -e 20200101 -m $MASKFILE"
opa_prex_or_die "cp sim_vs_clim_profiles/Fig_4.11*png Fig4.11"
opa_prex_or_die "cp sim_vs_clim_profiles/Fig_4.17*png Fig4.17"

DIR=static_clim
opa_prex_or_die "mkdir -p $DIR table4.6 table4.7 table4.10/ table4.11 table4.13/ table4.14"
# -------------------------------------------------------------------------

#python static_clim_validation_STD.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20170101 -e 20180101
#python static_clim_validation_STD_pH.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20170101 -e 20180101
opa_prex_or_die "python static_clim_validation.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20190101 -e 20200101"

DIR=static_clim_OpenSea
opa_prex_or_die "mkdir $DIR"
opa_prex_or_die "python static_clim_validation_OpenSea.py -i $STAT_PROFILES_DIR -o $DIR -m $MASKFILE -s 20190101 -e 20200101"

# -------------------------------------------------------------------------
opa_prex_or_die "cp $DIR/N1p-LAYER-Y-CLASS4-CLIM.txt $DIR/N3n-LAYER-Y-CLASS4-CLIM.txt table4.6/"
opa_prex_or_die "cp $DIR/O2o-LAYER-Y-CLASS4-CLIM.txt                                  table4.10/"

# PH-PROF-Y-CLASS4-CLIM-[BIAS/RMSD/CORR]-BASIN
# PCO-PROF-Y-CLASS4-CLIM-[BIAS/RMSD/CORR]-BASIN

#cp $DIR/Ac-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.13/
opa_prex_or_die "cp $DIR/ALK-LAYER-Y-CLASS4-CLIM.txt $DIR/DIC-LAYER-Y-CLASS4-CLIM.txt  table4.13/"
opa_prex_or_die "cp $DIR/pH-LAYER-Y-CLASS4-CLIM.txt  $DIR/pCO2-LAYER-Y-CLASS4-CLIM.txt table4.13/"

# ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14
# DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN calculated on 14 layers --> in table4.14

opa_prex_or_die "cp $DIR/N1p-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt $DIR/N3n-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.7"
opa_prex_or_die "cp $DIR/O2o-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt                                            table4.11"
#cp $DIR/Ac-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14
opa_prex_or_die "cp $DIR/ALK-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt  $DIR/DIC-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14"
opa_prex_or_die "cp $DIR/pH-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt   $DIR/pCO2-PROF-Y-CLASS4-CLIM-CORR-BASIN.txt table4.14"

### SECTION WITH NEW FIGURES FOR pCO2 COMPARISON/ANALYSIS ###
# Compare pCO2 vs SOCAT dataset: 
# Fig4.20  : PCO-SURF-M-CLASS4-CLIM-MEAN-BASIN
# table4.15: PCO-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt 
opa_prex_or_die "mkdir -p table4.15 Fig4.20 Fig4.21"
# First generate climatology tables:

MONTHLY_SURF=monthly_surf
opa_prex_or_die "mkdir -p $MONTHLY_SURF"
opa_prex_or_die "python monthly_clim_socat_pCO2.py"
opa_prex_or_die "python monthly_surf.py -i $STAT_PROFILES_DIR -o $MONTHLY_SURF -y 2019"
opa_prex_or_die "python table_pCO2vsSOCAT.py -i $MONTHLY_SURF"

opa_prex_or_die "mv pCO2-SURF-M-CLASS4-CLIM-RMSD-BASIN.txt table4.15/."
opa_prex_or_die "mv TOT_RMSD_pCO2vsSOCAT.txt table4.15/."

# generate new figure 4.20: comparison pCO2 BFM vs SOCAT dataset

opa_prex_or_die "python plot_month_pCO2vsSOCAT.py"
opa_prex_or_die "mv pCO2_monthly_tseries_Fig4.20.png Fig4.20/. "

# generate new figure 4.21: comparison pCO2 and CO2ariflux vs REAN (CLIM)
# Fig4.21: airseaCO2flux-SURF-M-CLASS4-CLIM-MEAN-BASIN

##REANDIR="/gpfs/scratch/userexternal/gbolzon0/REA/output/STAT_PROFILES/"
#REANDIR=/gpfs/scratch/userexternal/lfeudale/REANALYSIS/STAT_PROFILES/
#python plot_pCO2_CO2airflux_vs_REAN.py -r $REANDIR -i $STAT_PROFILES_MONTHLY_DIR -o Fig4.21/
#python plot_pCO2_CO2airflux_vs_REAN.py -r $REANDIR -i $STAT_PROFILES_DIR -o Fig4.21/

