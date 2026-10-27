#! /bin/bash
module load python/3.11.7--gcc--10.2.0

source $OPA_HOME/HOST/g100/venv/bin/activate
source ../../Sat/profile.inc

ONLINE_REPO=/g100_scratch/usera07ogs/a07ogs00/vsciortino/ONLINE/
OPA_ETCDIR=/g100_scratch/usera07ogs/a07ogs00/vsciortino/etc
TRAIN_DIR=$OPA_ETCDIR/PPCON_TRAIN
ONLINE_REPO_CLUSTERING=$ONLINE_REPO/PPCON/clustering/
mkdir - p $TRAIN_DIR
mkdir -p $ONLINE_REPO_CLUSTERING
echo "ONLINE_REPO: $ONLINE_REPO"
echo "OPA_ETCDIR: $OPA_ETCDIR"
echo "TRAIN_DIR: $TRAIN_DIR"
echo "ONLINE_REPO_CLUSTERING: $ONLINE_REPO_CLUSTERING"

# Set the local PPCON diff-list path inside the online repository.
LOCAL_INPUT=$ONLINE_REPO/PPCON/DIFF_floats.txt
LOCAL_INPUT_PY=$ONLINE_REPO/PPCON/DIFF_floats_py.txt
# Set the input update file path for the current run.
UPDATE_FILE=/g100_scratch/usera07ogs/a07ogs00/vsciortino/DIFF_floats.20261026-20:55:39.txt

###### PPCON generation based on UPDATE_FILE #############


> $LOCAL_INPUT
> $LOCAL_INPUT_PY

echo "Start selection of files using the shell scripting lines..."

while read -r line ; do
filename=$( echo $line | cut -d "," -f1);
  V1=${filename/coriolis\//}
  V2=${V1/profiles\//}
  superfloat_file=${ONLINE_REPO}/SUPERFLOAT/${V2}
  [ -f $superfloat_file ] && echo $line >> $LOCAL_INPUT

done < ${UPDATE_FILE}


echo "Start copying from SUPERFLOAT"
for filename in $( cat ${LOCAL_INPUT} | cut -d "," -f1 ) ; do
  V1=${filename/coriolis\//}
  V2=${V1/profiles\//}
   cp -v ${ONLINE_REPO}/SUPERFLOAT/${V2} ${ONLINE_REPO}/PPCON/${V2}
 done

echo "Start selection of files using the python script..."
my_prex_or_die "python -u reset_ppcon_dataset.py -i $UPDATE_FILE -o $LOCAL_INPUT_PY -I $ONLINE_REPO/SUPERFLOAT/ -O $ONLINE_REPO/PPCON/"

NR=$( wc -l $LOCAL_INPUT | awk '{print $1}')
echo "Number of lines in the local input file: $NR"
NR_PY=$( wc -l $LOCAL_INPUT_PY | awk '{print $1}')
echo "Number of lines in the local input file (Python): $NR_PY"

#if [[ $NR > 0 ]] ; then
#  my_prex_or_die "python -u clustering/clustering.py -i $ONLINE_REPO -u $LOCAL_INPUT -o $ONLINE_REPO_CLUSTERING"
#  my_prex_or_die "python -u make_generated_ds/generate_netcdf_netcdf4.py -t $ONLINE_REPO_CLUSTERING -m $TRAIN_DIR -p $ONLINE_REPO/PPCON"
#  my_prex_or_die "python dump_index.py -i $ONLINE_REPO/PPCON -o $ONLINE_REPO/PPCON/Float_Index.txt -t ppcon_float"
#fi