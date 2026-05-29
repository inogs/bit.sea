#! /bin/bash


. ../Sat/profile.inc

export ONLINE_REPO=/g100_scratch/usera07ogs/a07ogs00/V12C/ONLINE

DOWNLOAD_DIR=$ONLINE_REPO/CORIOLIS/download
LOCALDIR=$DOWNLOAD_DIR/tmp

UPDATE_FILE=DIFF_floats.$(date +\%Y\%m\%d-\%H:\%M:\%S).txt
UPDATE_FILE_PATH=$DOWNLOAD_DIR/$UPDATE_FILE
INDEXER_CORIOLIS=$DOWNLOAD_DIR/Med_floats.txt
INDEXER_CORIOLIS_OLD=$DOWNLOAD_DIR/Med_floats_OLD.txt

mkdir $LOCALDIR

#rename old syntetic_profile file and old output file
if [ -f "$INDEXER_CORIOLIS" ]; then
	my_prex "cp $INDEXER_CORIOLIS $INDEXER_CORIOLIS_OLD"
    echo "Old indexer file found, copying it to $INDEXER_CORIOLIS_OLD. Downloading only new profiles."
else
	my_prex "touch $INDEXER_CORIOLIS_OLD"
    echo "No old indexer file found, creating an empty one at $INDEXER_CORIOLIS_OLD. Downloading all profiles as new."
fi

#download new syntetic profile

REMOTE_PROFILE=/ifremer/argo/argo_synthetic-profile_index.txt.gz

my_prex_or_die "python $OPA_BITSEA/Float/argo_coriolis_updater.py -c $DOWNLOAD_DIR -R $REMOTE_PROFILE -O $INDEXER_CORIOLIS_OLD -u $UPDATE_FILE_PATH -i $INDEXER_CORIOLIS"
my_prex "cp $UPDATE_FILE_PATH $OPA_LOGDIR/daily/"

## download netcdf files in a tmp directory
CORIOLIS_DIR=$ONLINE_REPO/CORIOLIS
my_prex_or_die "python $OPA_BITSEA/Float/argo_downloader.py -c $CORIOLIS_DIR -u $UPDATE_FILE_PATH -i $INDEXER_CORIOLIS"

