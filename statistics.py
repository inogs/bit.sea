import numpy as np

def diff(Model, Ref):
    return Model - Ref

def number(Model):        
    return len(Model)

def variances(Model, Ref):
    return Ref.std()**2, Model.std()**2

def medians(Model, Ref):
    return np.median(Ref), np.median(Model)

def covariance(Model, Ref):
    array = (Model - Model.mean())*(Ref - Ref.mean())
    covariance = array.mean()
    return covariance

def covariance(Ref,Model,output_matrix=False):
    array = (Model - Model.mean())*(Ref - Ref.mean())
    return array.mean() if not output_matrix else array

def correlation(Model, Ref,output_matrix=False):
    array = (Model - Model.mean())*(Ref - Ref.mean())
    cov = array.mean() if not output_matrix else array
    return cov/(Model.std()*Ref.std())

def bias(Model, Ref):
    return (Model - Ref).mean()

def MSE(Model, Ref):
    '''Mean Square Error'''
    return ((Model - Ref)**2).mean()

def RMSE(Model, Ref):
    ''' Root mean Square Error'''
    return np.sqrt(((Model - Ref)**2).mean())
