def linear_shift(N_old,pres,shift,p_bot=600):
    '''
    Perform a linar shift of nitrate profile along the vertical direction starting "for default" from 600m depth upward
    Argument:
    * N_old * Original nitrate profile
    * shift * shift calculated with WOA at bottom (mean 600-1000m)
    * p_bot * depth from which starts the shift (for default is fixed to 600m)

    Return:
    * New_profile * Nitrate profile with a linear shift applied
    * qc          * qc modified to 8 (that means interpolated value)
    '''

    N_new=N_old-shift

    P600=pres[(pres>=p_bot)][0]
    New_profile = N_new.copy()
    Shift_Surf = N_old[0]-max(0.05,N_new[0])
    for iz, zz in enumerate(pres[(pres<=p_bot)]):
        New_profile[iz] = N_old[iz] - (Shift_Surf + (shift - Shift_Surf)*(pres[iz]-pres[0])/(P600-pres[0]))
        New_profile[iz]=max(0.05,New_profile[iz])  # Eliminate possible negative values


    return New_profile