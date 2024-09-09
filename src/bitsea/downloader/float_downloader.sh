#!/bin/bash

export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea
export ONLINE_REPO=/gpfs/scratch/userexternal/gbolzon0/V7C/ONLINE
DIR=$ONLINE_REPO/CORIOLIS/download
LOCALDIR=$ONLINE_REPO/CORIOLIS/download/tmp
DIFF_FILE=$DIR/DIFF_floats.txt
OPA_BITSEA=/gpfs/scratch/userexternal/gbolzon0/V7C/bit.sea
mkdir -p $DIR

#rename old syntetic_profile file and old output file

cp $DIR/argo_synthetic-profile_index_CORR.txt $DIR/argo_synthetic-profile_index_CORR_OLD.txt 

cp $DIR/Med_floats.txt $DIR/Med_floats_OLD.txt

#download new syntetic profile

HOST=ftp.ifremer.fr
USER=anonymous
REMOTEDIR=/ifremer/argo/

ncftpget -u $USER $HOST $DIR $REMOTEDIR/argo_synthetic-profile_index.txt.gz

gzip -dc $DIR/argo_synthetic-profile_index.txt.gz > $DIR/argo_synthetic-profile_index.txt

#correct file from blank parts

grep -v ",," $DIR/argo_synthetic-profile_index.txt | grep -v D.nc > $DIR/argo_synthetic-profile_index_CORR.txt

#launch the reader -> obtain the mediterranean floater file

python $OPA_BITSEA/Float/argo_reader.py -i $DIR/argo_synthetic-profile_index_CORR.txt -o $DIR/Med_floats.txt

#matching old vs new -> return a file called DIFF.txt in which there is the list of floats that are different between the two files

python $OPA_BITSEA/Float/argo_difference.py -N $DIR/Med_floats.txt -O $DIR/Med_floats_OLD.txt -o $DIFF_FILE



## download netcdf files in a tmp directory

REMOTEDIR=/ifremer/argo/dac/


mkdir -p $LOCALDIR
for I in `awk --field-separator=, '{print $1}'  $DIFF_FILE ` ; do
  ncftpget -u $USER $HOST $LOCALDIR $REMOTEDIR/$I
done


echo "Moving NetCDF files from TMP dir to ONLINE CORIOLIS file"
cd $LOCALDIR
for I in `ls *nc ` ; do
   WMO_DIR=${I:2:7}
   CORIOLIS_DIR=$ONLINE_REPO/CORIOLIS/${WMO_DIR}/
   mkdir -p $CORIOLIS_DIR
   mv $I $CORIOLIS_DIR
done
