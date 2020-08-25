import numpy as np
from datetime import datetime
from scipy import interpolate

gtime="20141209-08:45:00"

lat=17.6; lon=-24.3; pres=float(180); temp=16; psal=36.1; doxy=104
basedir="/Users/gbolzon/Documents/workspace/CANYON_B/CODES/CANYON_Training/"

d=datetime.strptime(gtime,"%Y%m%d-%H:%M:%S")
presgrid=np.loadtxt(basedir + 'CY_doy_pres_limit.csv',delimiter="\t")
doy=np.int(d.strftime('%j'))*360./365

if lon>180: lon=lon-360
x,y= np.meshgrid(presgrid[1:,0],presgrid[0,1:])
f = interpolate.interp2d(presgrid[1:,0],presgrid[0,1:],presgrid[1:,1:].T)
prespivot = f(lon,lat)

fsigmoid=1./(1 + np.exp((pres-prespivot)/50))
lonrad=lon*np.pi/180
doyrad=doy*np.pi/180

data=[lat/90, np.sin(lonrad), np.cos(lonrad), np.sin(doyrad)*fsigmoid, np.cos(doyrad)*fsigmoid, temp, psal, doxy, pres/2e4 + 1./((1+np.exp(-pres/300))**3)];
data = np.array(data)

poids=np.loadtxt(basedir + 'Fichier_poids_NO3_hidden_ascii20_17.sn')
poids=poids[:,2]


ne=9    # Number of inputs
nc1=20  # Number of neurons of the first hidden layer
nc2=17  # Number of neurons of the second hidden layer
ns=1    # Number of outputs


# WEIGHT AND BIAS PARAMETERS   
# weight and bias from the input layer to the first hidden layer w1, b1
b1=poids[0:nc1]
start=nc1+nc2+ns
end=start+ne*nc1
w1 = np.reshape(poids[start:end],(ne,nc1)).T

#% weight and bias from the first hidden layer to the second hidden layer
b2=poids[nc1:nc1+nc2]
start=nc1+nc2+ns+ne*nc1
end = start + nc1*nc2
w2=np.reshape(poids[start:end],(nc1,nc2)).T
#% weight and bias from the second hidden layer to the output layer w3, b3
start=nc1+nc2    
b3=poids[start:start+ns]
start=nc1+nc2+ns+ne*nc1+nc1*nc2
end = start + nc2*ns
w3=np.reshape(poids[start:end],(nc2,ns)).T

# Mean and standard deviation of the training dataset
# These values are used to normalize the inputs parameters
Moy  =np.loadtxt(basedir + 'moy_NO3.dat')
Ecart=np.loadtxt(basedir + 'std_NO3.dat')

# NORMALISATION OF THE INPUT PARAMETERS
rx=1 #[rx,~]=size(data);
data_N=(2./3)*(data-Moy[:ne])/Ecart[:ne]

# Two hidden layers
a=1.715905*np.tanh((2./3)*(np.dot(data_N,w1.T)+b1)) # input layer to first hidden layer
b=1.715905*np.tanh((2./3)*(np.dot(     a,w2.T)+b2)) # first hidden layer to second hidden layer
y=                         np.dot(     b,w3.T)+b3   # second hidden layer to output layer
# Y is the normalised output value of the neural network, 


# Denormalisation of the output of the NN for getting the true value 
y_rescaled=1.5*y*Ecart[ne]+Moy[ne]

# and put into same shape as the input variables
out=y_rescaled[0]
print out
