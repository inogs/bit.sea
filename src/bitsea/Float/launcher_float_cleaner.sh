#! /bin/bash
module load python/3.11.7--gcc--10.2.0

source $OPA_HOME/HOST/g100/venv/bin/activate

BITSEA_DIR=$OPA_HOME/HOST/g100/bit.sea/src/bitsea
source $BITSEA_DIR/Sat/profile.inc

ONLINE_REPO=/g100_scratch/usera07ogs/a07ogs00/vsciortino/ONLINE/

echo "ONLINE_REPO: $ONLINE_REPO"

FLOAT_DIR=$BITSEA_DIR/Float

###### SUPERFLOAT clean #############

my_prex "cd \"$ONLINE_REPO/SUPERFLOAT/\"; \ls -1R > SUPERFLOAT_contents.txt"

my_prex_or_die "python -u $FLOAT_DIR/float_cleaner.py -I $ONLINE_REPO/SUPERFLOAT/ -i $ONLINE_REPO/SUPERFLOAT/SUPERFLOAT_contents.txt"
my_prex_or_die "python -u $FLOAT_DIR/float_cleaner.py -I $ONLINE_REPO/SUPERFLOAT/ -i $ONLINE_REPO/SUPERFLOAT/SUPERFLOAT_contents.txt --force"

###### PPCON clean #############

my_prex "cd \"$ONLINE_REPO/PPCON/\"; \ls -1R > PPCON_contents.txt"

my_prex_or_die "python -u $FLOAT_DIR/float_cleaner.py -I $ONLINE_REPO/PPCON/ -i $ONLINE_REPO/PPCON/PPCON_contents.txt"
my_prex_or_die "python -u $FLOAT_DIR/float_cleaner.py -I $ONLINE_REPO/PPCON/ -i $ONLINE_REPO/PPCON/PPCON_contents.txt --force"
