#!/bin/bash

export PYTHONPATH=$PYTHONPATH:/gpfs/work/OGS20_PRACE_P/COPERNICUS/bit.sea

DIR=/galileo/home/userexternal/gcoidess/test_float_download/

#mkdir -p $DIR

#rename old syntetic_profile file and old output file

cp $DIR/argo_synthetic-profile_index_CORR.txt $DIR/argo_synthetic-profile_index_CORR_OLD.txt 

cp $DIR/Med_floats.txt $DIR/Med_floats_OLD.txt

#download new syntetic profile

HOST='ftp.ifremer.fr'
USER='anonymous'
REMOTEDIR=/ifremer/argo/

ncftpget -u $USER $HOST $DIR $REMOTEDIR/argo_synthetic-profile_index.txt.gz

gzip -dc $DIR/argo_synthetic-profile_index.txt.gz > $DIR/argo_synthetic-profile_index.txt

#correct file from blank parts

grep -v ",," $DIR/argo_synthetic-profile_index.txt | grep -v D.nc > $DIR/argo_synthetic-profile_index_CORR.txt

#launch the reader -> obtain the mediterranean floater file

python argo_reader.py -i $DIR/argo_synthetic-profile_index_CORR.txt -o $DIR/Med_floats.txt

#matching old vs new -> return a file called DIFF.txt in which there is the list of floats that are different between the two files

python argo_difference.py -N $DIR/Med_floats.txt -O $DIR/Med_floats_OLD.txt -o $DIR/DIFF_floats.txt

