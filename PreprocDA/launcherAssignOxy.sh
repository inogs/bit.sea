PVAR=15
INEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_Oxy/EOF_600_15sub_var$PVAR/
#INEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_N3n/EOF_600_15sub_var$PVAR/
#INEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_N3n/EOF_600_15sub_varzero/
OUTEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_Oxy/EOF_600_SUB_ALL_$PVAR/
#OUTEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_N3n/EOF_600_SUB_ALL_$PVAR/
#OUTEOF=/g100_scratch/userexternal/ateruzzi/DA_Oxy/EOF_N3n/EOF_600_SUB_ALL_varzero/
mkdir -p $OUTEOF

MASKMODEL=/g100_scratch/userexternal/camadio0/RUN_ChlNitOxy_SUM/wrkdir/MODEL/meshmask.nc
MASKVAR=/g100_scratch/userexternal/camadio0/RUN_ChlNitOxy_SUM/wrkdir/MODEL/DA_static_data/3D_VAR/GRID/BFM_gridFloat.nc
NLEVEOF=60

echo python assignEOFs_SUB_ALL.py -i $INEOF -o $OUTEOF -m $MASKMODEL -g $MASKVAR -n $NLEVEOF
