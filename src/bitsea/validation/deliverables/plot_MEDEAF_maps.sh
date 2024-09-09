#! /bin/bash
#SCR=/g100_work/OGS_prod100/OPA/V7C-prod/HOST/g100
SCR_HOME=/g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/INTERIM2021/bit.sea/

# mpirun -np 32 python /g100/home/usera07ogs/a07ogs00/OPA/V7C/HOST/g100/bitsea/build_layer_maps.py -l /g100/home/usera07ogs/a07ogs00/OPA/V7C/etc/static-data/POSTPROC/LogoEchoOGS4.png -o /g100_work/OGS_prod100/OPA/V7C/prod/wrkdir/forecast/1/POSTPROC/AVE_FREQ_1/MAP_GEN/MAP_ORIG -m /g100_work/OGS_prod100/OPA/V7C/prod/wrkdir/forecast/1/MODEL/meshmask.nc -i /g100_work/OGS_prod100/OPA/V7C/prod/wrkdir/forecast/1/MODEL/AVE_FREQ_1  -g ave*N1p.nc  -p /g100/home/usera07ogs/a07ogs00/OPA/V7C/HOST/g100/bitsea/postproc/Plotlist_bio.xml 

LOGO=/g100_work/OGS_prodC/OPA/V9C-prod/etc/static-data/POSTPROC/LogoEchoOGS4.png
PNG_DIR=/g100_scratch/userexternal/lfeudale/PNG_INTERIM/MAPS/
MASKFILE=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/MODEL/meshmask.nc
MASKFILE_CMCC=/gpfs/work/IscrC_REBIOMED/NRT_EAS6/PREPROC/MASK/ogstm/meshmask_CMCCfor_ogstm.nc
INPUTDIR=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/POSTPROC/output/MONTHLY_AVE/
INPUTDIR_PHYS=/g100_scratch/userexternal/gcoidess/REANALISI_INTERIM/wrkdir/POSTPROC/output/MONTHLY_PHYS/

#PLOTLIST=/g100_work/OGS_prod100/OPA/V7C-prod/HOST/g100/bitsea/postproc/Plotlist_bio.xml
PLOTLIST=/g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/INTERIM2021/bit.sea/postproc/Plotlist_bio_reduced.xml 
PLOTLIST_PHYS=/g100_work/OGS_prod100/OPA/V7C-prod/HOST/g100/bitsea/postproc/Plotlist_phys.xml 

cd $SCR_HOME

# NOTE: set the correct time interval in "build_layer_maps_INTERVAL.py"
# MAPS FOR BGC VARS:
python /g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/INTERIM2021/bit.sea/build_layer_maps_INTERVAL.py -l $LOGO -o $PNG_DIR -m $MASKFILE -i $INPUTDIR -g ave*N1p.nc -p $PLOTLIST

# MAPS FOR PHYS VARS:
python /g100_scratch/userexternal/lfeudale/REANALYSIS/REAN24/validation/INTERIM2021/bit.sea/build_layer_maps_INTERVAL.py -l $LOGO -o $PNG_DIR -m $MASKFILE -i $INPUTDIR_PHYS -g ave*votemper.nc -p $PLOTLIST_PHYS

cd $SCR_HOME/validation/deliverables
