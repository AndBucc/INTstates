# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:18:56 2023

@author: bucca
"""
# Cramer's V (Autocorrelation Function for categorical time series data)

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder

# Load your categorical series here. 
#%% load data
import scipy.io as sio
maptimeseries = group_data #here, the output of MATLAB's INT state map time series
all_cramers = []

#%% rest of code
for dat in range(mapsarray.shape[0]):
        
    data_beforemakeup = mapsarray[dat,:]-1
    label_encoder = LabelEncoder()
    data = label_encoder.fit_transform(data_beforemakeup)
    
    Tlen = len(data)
    
   
    states = np.unique(data)
    nostates = len(states)
    
    #nostates = 7 #if you want to input manually, uncomment
    #datanum = np.array([np.where(states == s)[0][0] for s in data])
    
    # Binarization
    bincodes = np.eye(nostates)
    databin = bincodes[data] #when input is already numerical
    #databin = bincodes[datanum] #when input is A,B,C, etc...
    
    # Relative frequencies
    hatpi = np.mean(databin, axis=0)
    
    maxlag = 100
    hatbivprob = np.zeros((nostates, nostates, maxlag))
    
    for k in range(1, maxlag + 1):  # For each lag
        for i in range(nostates):  # For each lagged vector representing a category
            for j in range(nostates):  # For each vector representing a category
                hatbivprob[i, j, k-1] = np.mean(databin[k:Tlen, i] * databin[0:Tlen-k, j])
    
    # Compare
    indprob = np.outer(hatpi, hatpi)
    
    # Cramer's V
    cramer = np.zeros(maxlag)
    for k in range(maxlag):
        cramer[k] = np.sqrt(np.sum((hatbivprob[:, :, k] - indprob)**2 / indprob) / (nostates - 1))
    all_cramers.append(cramer)
    #peak_cramer[dat] = 
    
    # Plot Cramer's V
    # plt.plot(range(1, maxlag + 1), cramer, marker='o', linestyle='-', linewidth=2)
    # plt.xlabel('k (Lag)')
    # plt.ylabel("Cramer's V(k)")
    # plt.title("Cramer's V for Different Lags")
    # plt.ylim(-1, 1)
    # plt.grid(True)
    # plt.show()

#%% IGNORE IF YOU HAVE ALREADY IMPORTED "decayrate.py"
# fit an exponential decay to Cramer's results (the whole pipeline)

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Exponential decay function for curve fitting
# y = A e^(-dr*t) + C. 
# x: x-axis, dr: decay rate (>0), (A + C): value of y when t = 0
def expdecay(x, dr, A, C):
    return A * np.exp(-dr * x) + C


def decayrate(xdata, ydata, trim=True, trimtp=5, dr_guess = 0.5, A_guess = 0.5, C_guess = 0.25):
    """
    Parameters
    ----------
    xdata : 1D numpy array of your x axis
    ydata : 1D numpy array of your y axis
    
    trim : Logical
        If True, trim what comes after the y axis reaches 0.
    trimtp : Integer
        Number of time points to count after ACF reaches 0 before trimming.
        For example, if this is 5, the data will be trimmed after 5 lags from touching 0
        point.
    
    *_guess: Initial guess (starting point) in nonlinear optimization.
    Note that this is a bit arbitrary. What matters is that this has to be
    same in all subjects. I provided some guesses as initial points, I remember that your 
    plots were starting from 0.75 and ending aroung 0.25. If you set C = 0.25 (asymptotical ending), 
    A = 0.50 when t = 0. If t doesn't start from 0, you can estimate initial guess of A by solving 
    the above equation for A when t = your starting point on x axis and c = 0.25 (since it seems to be 
    the asymptotical value).

    Note that exponential decay theoretically never reaches 0. I used the trim to make the estimations
    more stable. Play around with it, if it doesn't work, set it to false.

    Returns
    -------
    dr : decay rate of estimated exponential function

    """

    # Decoration: Get rid of the portion after ACW-0.
    lags = []
    if trim:
        acw_0 = np.argmax(ydata <= 0.25)
        ydata = ydata[0 : (acw_0 + trimtp)]
        lags = lags[0 : (acw_0 + trimtp)]

    params, _ = curve_fit(expdecay, xdata, ydata, p0=(dr_guess, A_guess, C_guess), bounds=(0.0, np.inf), maxfev=5000)
    dr = params[0]
    A = params[1]
    C = params[2]
    return dr, A, C


# Quality control function
def check_fit(dr, A, C, xdata, ydata):
    """ Plot to check the goodness of fit"""
    plt.plot(xdata, expdecay(xdata, dr, A, C))
    plt.plot(xdata, ydata)
    plt.show()


#%% actual fitting. Example with a group of 10 subjects. Substitute with your data.
decays_sample = []
for i in range(len(all_cramers)):
    dr = decayrate(np.arange(1,10,1),all_cramers[i][0:30], trim = False)
    decays_sample.append(dr)
    

for i in range(10):
    temp = decays_sample[i]
    check_fit(temp[0],temp[1],temp[2],np.arange(1,10,1), all_cramers[i][0:30])
    
    
d_sample = [decays_sample[i][0] for i in range(10)]
    
