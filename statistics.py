import numpy as np

def diff(Model, Ref):
    return Model - Ref

def number(Model):        
    return len(Model)

def variances(Model, Ref):
    # Check numpy variance function so as not
    # to do the **2 operation
    return np.nanvar(Ref),np.nanvar(Model)
    #return Ref.std()**2, Model.std()**2

def medians(Model, Ref):
    return np.nanmedian(Ref), np.nanmedian(Model)

def covariance(Model,Ref,output_matrix=False):
    # Means of Ref and Model taking into account that
    # there might be NaN's inside
    Model_mean = np.nanmean(Model)
    Ref_mean   = np.nanmean(Ref)
    cov_array  = (Model - Model_mean)*(Ref - Ref_mean)
    return np.nanmean(cov_array) if not output_matrix else cov_array

def correlation(Model,Ref,output_matrix=False):
    # Use the covariance function to obtain the covariance
    cov = covariance(Model,Ref,output_matrix)
    Model_std = np.nanstd(Model)
    Ref_std   = np.nanstd(Ref)
    return cov/(Model_std*Ref_std)

def bias(Model, Ref):
    return np.nanmean(Model - Ref)

def MSE(Model, Ref):
    '''Mean Square Error'''
    return np.nanmean((Model - Ref)*(Model - Ref))

def RMSE(Model, Ref):
    ''' Root mean Square Error'''
    return np.sqrt(MSE(Model,Ref))
