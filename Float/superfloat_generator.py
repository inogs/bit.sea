from commons import calculated_depths
import numpy as np
import scipy.io.netcdf as NC
import os


def general_quenching(profile_obj, Pres, Profile, Qc):
    '''
    Manages the quenching for Coriolis profiles
    When mld > 80m (flag 1 ) it doesn't do anything
    When Z_top > mld + 5m (flag 2 ) idem
    In the remaining cases: it applies the quenching, from mld to top

    In case of cutted profile (Z_top > 10.) the quench value is replied up to the sea surface
    on a set of extrapolation points, every 5 meters.
    '''
    PresT, Temp, _ = profile_obj.read('TEMP', read_adjusted=False)
    mld = calculated_depths.mld(Temp, PresT, zref=0, deltaTemp=0.1)
    if mld < 80:
        if Pres[0] <= mld+5:
            ii = Pres <= mld + 5
            descending_ordered = np.sort(Profile[ii])[-1::-1]
            n_chosen = min(3, len(descending_ordered))
            if n_chosen==0: n_chosen = 1
            quench_value = descending_ordered[:n_chosen].mean()
            ii = Pres < mld
            Profile[ii] = quench_value

            if Pres[0] > 10. : # cutted profiles
                creation_step = 5 #m
                Top_Pressures= np.arange(0.5,Pres[0],creation_step)
                Top_Values   = np.ones_like(Top_Pressures)* quench_value
                Top_QC       = np.ones_like(Top_Pressures)* 1

                Pres    = np.concatenate((Top_Pressures, Pres), axis=0)
                Profile = np.concatenate((Top_Values, Profile), axis=0)
                Qc      = np.concatenate((Top_QC,          Qc), axis=0)

    return Pres, Profile, Qc


def quenching(profile_obj, PresChl, Chl, chl_lov_zero):
    '''
    Apply quenching by substitution of chl values
    chl(Mixed Layer) = max(chl(Mixed Layer))


    Arguments:
    * profile_obj  * a Coriolis profile object
    * Pres_chl     * numpy 1D-array of coriolis pressure
    * Chl          * numpy 1D-array of coriolis chl
    * chl_lov_zero * value of LOV quench; an algorithm will try to guess the Mixed Layer Depth
                     calculated by LOV, in order to use it here to calculate quench value
                     if None, the Coriols Mixed Layer Depth will be used.
    Returns:
    * CHL          * numpy 1D-array of corrected chlorophyll
    '''
    Quenched_Chl = Chl.copy()
    if profile_obj.time.month in range(1,13):
        PresT, Temp, _ = profile_obj.read('TEMP', read_adjusted=False)
        mld = calculated_depths.mld(Temp, PresT, zref=0, deltaTemp=0.1)
        if mld > 200:mld=20

        ii=PresChl<=mld
        base_mld_level=ii.sum()
        if chl_lov_zero is not None:
            nTries = 25
            if  nTries > len(PresChl)-base_mld_level-1:
                nTries = len(PresChl)-base_mld_level-1

            quench_value = np.zeros((nTries,), np.float32)

            diff_fromlov = np.zeros((nTries,), np.float32)
            for down_level in range(nTries):
                nLev = base_mld_level + down_level
                if nLev==0: nLev=1
                quench_value[down_level] = Chl[:nLev].max()
                diff_fromlov[down_level] = np.abs(quench_value[down_level]-chl_lov_zero)

            k = diff_fromlov.argmin()
            nLev= base_mld_level + k
            Quenched_Chl[:nLev] = quench_value[k]
        else:
            nLev = base_mld_level + 1
            quench_value = Chl[:nLev].max()
            Quenched_Chl[:nLev] = quench_value
    return Quenched_Chl

def are_identical(chl1, chl2, threshold=0.05):
    diff = chl1 - chl2
    if np.abs(diff).max()< threshold:
        return True
    else:
        return False
def are_shifted(chl1, chl2):
    diff = chl1 - chl2
    l = np.median(diff) - np.mean(diff)
    if l<1.e-3:
        return True
    else:
        return False
def has_drift(Pres,Profile):
    if Pres.max()<200:
        return False
    i200 = np.abs(Pres-200).argmin()
    Delta = Profile[-1] - Profile[i200]
    if np.abs(Delta) > 0.05:
        return True
    else:
        return False

def synthesis_profile(pLov, pCor):
    '''
    Arguments:
    * pLov, pCor * two profile objects
    Returns
    * Pres, Value, Qc * numpy ndarrays
    '''

    PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
    if ('CHLA' in pCor.available_params) and  (pCor._my_float.status_var('CHLA')=='A'):
        PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=True)
    else:
        return PresL, ValueL, QcL
    if len(PresC)<5:

        if len(PresL)>5:
            print "few values in Coriolis: using LOV"
            return PresL, ValueL, QcL
        else:
            print "few values in either dataset for " + pLov._my_float.filename
            return None, None, None
    if has_drift(PresC, ValueC):
        print pCor._my_float.filename + " has drift"
        return None, None, None

    if len(ValueL)> 0:
        chl_lov_zero=ValueL[0]
    else:
        chl_lov_zero = None
    CHL = quenching(pCor, PresC, ValueC,chl_lov_zero)

    ii=(PresC >= 400) & (PresC <= 600)
    if ii.sum() > 0:
        shift = CHL[ii].mean()
        CHL = CHL - shift
    return PresC, CHL, QcC



def exist_valid(filename):
    '''Returns true if file exists and is a valid NetCDF file'''
    if not os.path.exists(filename):
        return False
    good = True
    try:
        ncIN = NC.netcdf_file(filename,'r')
    except:
        good = False
    if good:
        ncIN.close()
        return True
    else:
        return False

def exist_variable(variable, filename):
    ncIN = NC.netcdf_file(filename,'r')
    variables=ncIN.variables.keys()
    ncIN.close()
    return variable in variables

def exist_valid_variable(variable, filename):
    if not exist_valid(filename): return False
    return exist_variable(variable, filename)

def writing_mode(filename):
    '''
    Returns 'a' for valid NetCDF files, else 'w' (overwrite)
    '''
    writingmode='w'
    if exist_valid(filename): writingmode='a'
    return writingmode

class Metadata():
    def __init__(self, filename):
        self.filename = filename
        self.status_var = 'n'


def read_float_update(input_file):
    mydtype= np.dtype([
        ('file_name','S200'),
        ('date','S200'),
        ('latitude',np.float32),
        ('longitude',np.float32),
        ('ocean','S10'),
        ('profiler_type',np.int),
        ('institution','S10'),
        ('parameters','S200'),
        ('parameter_data_mode','S100'),
        ('date_update','S200')] )

    INDEX_FILE=np.loadtxt(input_file,dtype=mydtype, delimiter=",",ndmin=1,skiprows=0)

    return INDEX_FILE



if __name__=="__main__":
    from instruments import bio_float
    from instruments import lovbio_float
    from commons.time_interval import TimeInterval
    from basins.region import Rectangle
    import matplotlib.pyplot as pl
    import sys
    TI = TimeInterval('2012','2020','%Y')
    R = Rectangle(-6,36,30,46)

#     PROFILES_LOV =lovbio_float.FloatSelector('BBP700', TI, R)
#     for ip, pLov in enumerate(PROFILES_LOV[:]):
#         pCor=bio_float.from_lov_profile(pLov)
#         if pCor is not None:
#             PresL, ValueL, QcL = pLov.read('BBP700',read_adjusted=False)
#             if 'BBP700' in pCor.available_params:
#                 PresC, ValueC, QcC = pCor.read('BBP700',read_adjusted=False)
#                 if len(ValueC)==len(ValueL):
#                     if not are_identical(ValueC, ValueL, threshold=1.e-3):
#                         print "different", ip
#                         fig,ax =pl.subplots()
#                         ax.plot(ValueL,PresL,'r', label="LOV")
#                         ax.plot(ValueC  , PresC,'b.-', label="COR")
#                         ax.invert_yaxis()
#                         ax.grid()
#                         ax.legend()
#                         fig.show()
#                         sys.exit()
#                         if are_shifted(ValueC, ValueL):
#                             pass
#                         sys.exit()
#                 else:
#                     print "different lengths", ip, len(ValueC), len(ValueL)
#
#
#     sys.exit()
#
#
#     PROFILES_LOV =lovbio_float.FloatSelector('PAR', TI, R)
#     for ip, pLov in enumerate(PROFILES_LOV[:1000]):
#         pCor=bio_float.from_lov_profile(pLov)
#         if pCor is not None:
#             PresL, ValueL, QcL = pLov.read('PAR',read_adjusted=False)
#             if 'DOWNWELLING_PAR' in pCor.available_params:
#                 PresC, ValueC, QcC = pCor.read('DOWNWELLING_PAR',read_adjusted=False)
#                 if not are_identical(ValueC, ValueL, 10):
#                     print "different", ip
#                     fig,ax =pl.subplots()
#                     ax.plot(ValueL,PresL,'r', label="LOV")
#                     ax.plot(ValueC  , PresC,'b.-', label="COR")
#                     ax.invert_yaxis()
#                     ax.grid()
#                     ax.legend()
#                     fig.show()
#                     # sys.exit()
#
#     sys.exit()
    PROFILES_LOV =lovbio_float.FloatSelector('CHLA', TI, R)   
#    PROFILES_LOV = lovbio_float.filter_by_wmo(PROFILES_LOV, "6901769")
    counter=0
    BAD_LIST=[]
    for ip, pLov in enumerate(PROFILES_LOV[:]):
        #if pLov._my_float.wmo == "6901765" : continue # save time
        #if pLov._my_float.wmo != "6901772" : continue # save time
        pCor=bio_float.from_lov_profile(pLov)
        
        if pCor is not None:
            PresL, ValueL, QcL = pLov.read('CHLA',read_adjusted=True)
            if 'CHLA' in pCor.available_params:
                if pCor._my_float.status_var('CHLA')=='A':
                    PresC, ValueC, QcC = pCor.read('CHLA',read_adjusted=True)
                else:
                    print pCor._my_float.filename + " without CHLA_ADJUSTED"
                    continue
            else:
                print pCor._my_float.filename + " without CHLA"
                continue

            if len(ValueC) < 5:
                print "few values in Coriolis for " + pCor._my_float.filename
                continue

            if len(ValueL)<  1:
                print "no LOV data"
                continue

            if has_drift(PresC, ValueC):
                print pCor._my_float.filename + " has drift"
                continue

            ValueCorig = ValueC.copy()
            PresCorig  = PresC.copy()
            PresC, CHL, Qc = general_quenching(pCor, PresC, ValueC, QcC)


            ii=(PresC >= 400) & (PresC <= 600)
            if ii.sum() > 0:
                shift = CHL[ii].mean()
                CHL = CHL - shift


            if True: #not are_identical(CHL, ValueL):
                if True: # not are_shifted(CHL, ValueL):
                    BAD_LIST.append(ip)
                    fig,ax =pl.subplots()
                    ax.plot(ValueL,PresL,'r', label="LOV")
                    ax.plot(ValueCorig  , PresCorig,'m.-', label="CORorig")
                    ax.plot(CHL  , PresC,'b.-', label="COR")
                    ax.invert_yaxis()
                    ax.grid()
                    ax.legend()
                    fig.show()
                    counter+=1
                    #pl.savefig(pCor.ID()+ ".png")
                    #pl.close(fig)
                    if counter > 10: sys.exit()
