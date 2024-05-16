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
        acw_0 = np.argmax(ydata <= 0.0)
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

