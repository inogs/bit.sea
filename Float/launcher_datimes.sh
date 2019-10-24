
DAFREQ=1
#VARDA=N3n
#VARDA=P_l
DA_HH=12


for VARDA in N3n P_l; do
echo python  profile_dates_DAfreq.py -f $DAFREQ -v $VARDA -t $DA_HH
echo ln -s daTimes_floatfreq${DAFREQ}_${VARDA} daTimes_float_${VARDA}
done

# Modify varlist in merge_daTimes.py
echo python merge_daTimes.py
