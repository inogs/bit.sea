# code to reconstruct profiles in NRT:

#cloned by the git of GloriaPP branch "CODE_loss_attention_max_PPCon" renamend by Amadio "CODE_loss_attention_max_PPCon_CA"

# CODE_loss_attention_max_PPCon --> it is the NNmodel  branch which reconstruct profiles 
# in CODE_loss_attention_max_PPCon/results  the trained NN model outputs

# the CODE works when:         in      ONLINE REPO/ 
                         I have a                 /SUPERFLOAT 
                         and it create            /PPCON 


# the CODE's input is: CASE  A)     Float_index        --> all the files are written
                       CASE  B)     Diff_floats[$DATE] --> a specific date is updated  
 

the CODE main programs are: 1) clustering.py         --> it create a tensor of all info: ds_sf_clustering.csv
 				   |
				   |   CASE A)  cp SUPERFLOAT/**/*nc in PPCON 
				   |   CASE B)  cp SUPERFLOAT[SUBSET_DATE] in PPCON
 				   V
                            2) make_generated_ds/generate_netcdf_netcdf4.py   --> it append the NN-Rec vars in a PPCON dataset


# 3 launcher to have:
Launcher_allsuperfloat_only.sh ---> entire dataset PPCON from SUPERFLOAT
Launcher_NRT_updt_ONLY.sh      ---> update PPCON from SUPERFLOAT updated NRT
Launcher_validation_NRT.sh     ---> to launch validation and create profiles as in medeaf website
