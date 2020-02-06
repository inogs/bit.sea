
from instruments import optbio_float_2019
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from statistics import *


def euphotic(Pres, Ed):
    EUPH_model         = [x for x in range(len(Pres)) if Ed.Model[x]/Ed.Model[0] > 0.37]
    EUPH_float         = [x for x in range(len(Pres)) if Ed.Ref[x]/Ed.Ref[0] > 0.37]
    ind_model          = EUPH_model[-1]
    ind_float          = EUPH_float[-1]
    ind  = ind_float  if ind_float < ind_model else ind_model
    EUPH = EUPH_float if ind_float < ind_model else EUPH_model
    zeu                = Pres[ind]
    Pres_Eu            = Pres[EUPH]
    Ed_Eu_MODEL        = Ed.Model[EUPH]
    Ed_Eu_FLOAT        = Ed.Ref[EUPH]
    return zeu, Pres_Eu, Ed_Eu_MODEL, Ed_Eu_FLOAT


def calc_KD(varname, requestor, basin, TI, M, mydepth, var_arg):

	Profilelist=optbio_float_2019.FloatSelector(varname,requestor, basin) # This is where you define the time requestor and the subbasin

	MODEL = []
	FLOAT = []

	Kd_Model = 0.1
	Kd_Float = 0.1

	for p in Profilelist:
	    profile_ID = p.ID()
	    
	    #phase 1. Read BGC-ARGO profiles
	    Pres, Ed_float,  Qc = p.read(varname)  
	    Lon = p.lon
	    Lat = p.lat
	    timestr = p.time.strftime("%Y%m%d-%H:%M:%S")
	    nLevels = len(Pres)

	    if not TI.contains(p.time):
	        continue
	    
	    #phase 2. Read model data - in theory you can use this also for float data
	    Ed_matchup = M.getMatchups([p], mydepth, var_arg)

	    #phase 3. Euphotic depth range calculation
	    zmax, PresEu, Ed_MODEL, Ed_REF = euphotic(Pres, Ed_matchup)

	    
	    #phase 4. Calculate Kd
	    if len(PresEu) < 5: continue
	    if PresEu[1] - PresEu[0] > 9. : continue
	    if PresEu[4] > 15. : continue

	    func = lambda z, Kd: np.exp(-Kd*z)

	    poptM, pcovM = curve_fit(func, PresEu-PresEu[0], Ed_MODEL/Ed_MODEL[0], p0=Kd_Model)#,bounds=(0.001,0.005))
	    poptF, pcovF = curve_fit(func, PresEu-PresEu[0], Ed_REF/Ed_REF[0],     p0=Kd_Float)#,bounds=(0.001,0.005))

	    Kd_Model = poptM[0]   ; Kd_Float = poptF[0]

	    MODEL.append(Kd_Model)
	    FLOAT.append(Kd_Float)

	    if Kd_Model < 0.01:
        	continue

	return np.array(FLOAT), np.array(MODEL)

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


