#! /bin/bash


. ../Sat/profile.inc

export ONLINE_REPO=/g100_scratch/usera07ogs/a07ogs00/V12C/ONLINE
CORIOLIS_DIR=$ONLINE_REPO/CORIOLIS
#REMOTEDIR=/ifremer/argo/dac/ # no more utilized argument
UPDATE_FILE=/g100_work/OGS_prod2528/OPA/V12C-prod/log/daily/DIFF_floats.20260424-20:53:57.txt
INDEXER_CORIOLIS=$ONLINE_REPO/CORIOLIS/download/Med_floats.txt
#NCFTPGET_PATH=/leonardo/home/usera07ogs/a07ogs00/OPA/V12C-dev/HOST/leonardo/bin/ncftpget

# version that uses ncftpget
#my_prex_or_die "python $OPA_BITSEA/Float/argo_downloader.py -c $CORIOLIS_DIR -u $UPDATE_FILE -i $INDEXER_CORIOLIS" # optional argument: -g $NCFTPGET_PATH

my_prex_or_die "python $OPA_BITSEA/Float/argo_downloader_2.py -c $CORIOLIS_DIR -u $UPDATE_FILE -i $INDEXER_CORIOLIS"

