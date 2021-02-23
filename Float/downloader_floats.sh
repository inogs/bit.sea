#!/bin/bash

DIR=/galileo/home/userexternal/gcoidess/test_float_download/
HOST='ftp.ifremer.fr'
USER='anonymous'
REMOTEDIR=/ifremer/argo/dac/
LOCALDIR=$DIR/download/
mkdir -p $LOCALDIR

for I in `awk --field-separator=, '{print $1}' $DIR/DIFF_floats.txt`
do ncftpget -u $USER $HOST $LOCALDIR $REMOTEDIR/$I
done

