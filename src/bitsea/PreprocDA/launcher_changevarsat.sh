

MESHMASK=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/meshmask_v12c.nc
INDIR=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/VARSAT_V9C/ # Copiati da /g100_work/OGS_devC/V9C/RUNS_SETUP/PREPROC/DA/DA_static_data (sono quelli usati nella QUID_NECCTON)
# LIMERRCOAST=0.1 #Valore limite dell'errore usato in costa
LIMERRCOAST=nan # Mettendo nan NON corregge in costa (testato)
FIRSMONTH=6 # primo mese d'estate
LASTMONT=8 # ultimo mese d'estate
ADDERRSUMMER=0.005 # Mettendo nana NON dovrebbe correggere in costa (NON TESTATO)
SATMAP=Mesh24

OUTDIR=/g100_scratch/userexternal/ateruzzi/PerCarolina/ADJUST_VARSAT/VARSATadj_COAST${LIMERRCOAST}_SUMMER${ADDERRSUMMER}

mkdir -p $OUTDIR

echo python adjust_varsat -i $INDIR -o $OUTDIR -m $MESHMASK -s $FIRSMONTH -e $LASTMONT -v $ADDERRSUMMER -c $LIMERRCOAST -mm $SATMAP
