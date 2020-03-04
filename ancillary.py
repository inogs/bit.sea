
from instruments import optbio_float_2019
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from statistics import *


def find_DCM(CHL, z):
#	# Filter low values
	CHL_filt = CHL[CHL>0.1]
	z_filt   = z[CHL>0.1]
	if (len(CHL_filt) == 0): return np.nan, np.nan
	# Find maximum
	DCM_ind = np.argmax(CHL_filt)
	DCM_z   = z_filt[  DCM_ind]
	DCM_val = CHL_filt[DCM_ind]
	# Match conditions
	if DCM_z < 30. or DCM_z > 200.: return np.nan, np.nan
	#if DCM_val/CHL[0]< 1.5:      return np.nan, np.nan
	return DCM_z, DCM_val

def calc_DCM(varname, requestor, basin, TI, M, mydepth, var_arg):

	Profilelist=optbio_float_2019.FloatSelector(varname,requestor, basin) # This is where you define the time requestor and the subbasin

	MODEL_DCM = []    ;   MODEL_CHL_DCM = []
	FLOAT_DCM = []    ;   FLOAT_CHL_DCM = []

	for p in Profilelist:
	    profile_ID = p.ID()
	    
	    #phase 1. Read BGC-ARGO profiles
	    Pres, Chl_float,  Qc = p.read(varname)  
	    Lon = p.lon
	    Lat = p.lat
	    timestr = p.time.strftime("%Y%m%d-%H:%M:%S")
	    nLevels = len(Pres)

	    if not TI.contains(p.time):
	        continue
	    
	    #phase 2. Read model data - in theory you can use this also for float data
	    Chl_matchup = M.getMatchups([p], mydepth, var_arg)

	    #phase 3. Calculate DCM depth and val

	    DCM_float, CHL_DCM_float = find_DCM(Chl_float, Pres) # Or Chl_matchup.Ref
	    DCM_model, CHL_DCM_model = find_DCM(Chl_matchup.Model, Pres)

	    MODEL_DCM.append(DCM_model)
	    FLOAT_DCM.append(DCM_float)

	    MODEL_CHL_DCM.append(CHL_DCM_model)
	    FLOAT_CHL_DCM.append(CHL_DCM_float)

	return np.array(FLOAT_DCM), np.array(MODEL_DCM), np.array(FLOAT_CHL_DCM), np.array(MODEL_CHL_DCM)


def calc_statistics(FLOAT, MODEL):

	count      = number(MODEL)
	corr_coeff = correlation(MODEL, FLOAT,output_matrix=False)
	bias_val   = bias(MODEL, FLOAT)
	slope, intercept, r_value, p_value, std_err = stats.linregress(FLOAT,MODEL)
	sigma      = RMSE(MODEL, FLOAT)
	a          = intercept
	b          = slope
	float_mean = FLOAT.mean()
	model_mean = MODEL.mean()
	return count, sigma, bias_val, corr_coeff, r_value, p_value, float_mean, model_mean  


def plot_basin(iseas, axs, FLOAT, MODEL, color, basin, titlestr, xpos, ypos, sigma, bias_val, corr_coeff,b,a,count):

    if iseas in range(2):
        index_row = 0
        index_col = iseas
    else:
        index_row = 1
        index_col = iseas-2
        
    axs[index_row, index_col].scatter(FLOAT, MODEL, marker='o', s=5.0, c=color,  label=basin.extended_name[0:4])
    axs[index_row, index_col].set_xlabel('Kd float [$m^{-1}$]')
    axs[index_row, index_col].set_ylabel('Kd model [$m^{-1}$]')
    axs[index_row, index_col].legend(loc='best')

    x_max      = np.max(FLOAT) *  1.1
    x_reg      = np.arange(0., x_max)
    axs[index_row, index_col].plot(x_reg,a+b*x_reg,color)
    axs[index_row, index_col].plot(x_reg,x_reg,'k--')
    
    textstr='$\mathrm{RMS}=%.2f$\n$\mathrm{Bias}=%.2f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2f$'%(sigma, bias_val, corr_coeff,b,a,count)
    
    axs[index_row, index_col].text(xpos, ypos, textstr, transform=axs[index_row, index_col].transAxes, color = color, fontsize=8, verticalalignment='top',bbox=dict(facecolor='white', alpha = 0.5, edgecolor=color))
    axs[index_row, index_col].set_title(titlestr[iseas])
    
    axs[index_row, index_col].legend(loc='lower center', bbox_to_anchor=(0.5, 0.85), ncol=2, fancybox=False, shadow=False).get_frame().set_alpha(0.5)
    
    return axs

def writefile(filestat, var_arg, count, sigma, corr_coeff, bias_val, r_value, p_value, float_mean, model_mean):
    fid = open(filestat,'wb')
    
    fid.write("%s %.2f %.2f %.2f %.2f %.2f %.2f %.2f  \n" % (var_arg, sigma, bias_val, corr_coeff, r_value, p_value, float_mean, model_mean) )

    fid.close()
    return fid


