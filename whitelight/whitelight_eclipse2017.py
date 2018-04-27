import sys

import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange
from tqdm import tqdm
import pandas as pd

import occultnl
from wave_solution import orbits


def lightcurve(p, x, y, err, sh, rprs,transit=False,fjac=None):
    """ Function used by MPFIT to fit data to lightcurve model. 
    
    Inputs: p: input parameters that we are fitting for
    x: Date of each observation, to be converted to phase
    y: Flux of each time series observation
    err: Error on flux of each time series point
    sh: Parameter for wavelength shift on the detector (for each exposure)
    rprs: Ratio of planet to star radii
    trasnit: True for transit, false for eclipse

    Output: Returns weighted deviations between model and data to be minimized
    by MPFIT. """
    
    
    # params= [fp,flux0,epoch,m,HSTP1,HSTP2,HSTP3,HSTP4,HSTP5,HSTP6,xshift1
    # ,xshift2,xshift3,xshift4,xshift5,xshift6,inclin,MsMpR,c1,c2,c3,c4,Per,T0]
    Per = p[22]  
    JD = 2400000.5             
    Gr = 6.67259e-11
    HSTper = 96.36 / (24.*60.)    
    inclin = p[16]
    MsMpR = p[17]

    phase = (x-p[2])/Per 
    phase = phase - np.floor(phase)
    phase[phase > 0.5] = phase[phase > 0.5] -1.0

    HSTphase = (x-p[23])/HSTper
    HSTphase = HSTphase - np.floor(HSTphase)
    HSTphase[HSTphase > 0.5] = HSTphase[HSTphase > 0.5] -1.0


    systematic_model_mpfit = ((phase*p[3] + 1.0) 
                              * (HSTphase*p[4] + HSTphase**2.*p[5] + HSTphase**3.*p[6]
                                 + HSTphase**4.*p[7] + HSTphase**5.*p[8] + HSTphase**6.*p[9] + 1.0)
                              * (sh*p[10] + sh**2.*p[11] + sh**3.*p[12] + sh**4.*p[13]
                                 + sh**5.*p[14] + sh**6.*p[15] + 1.0))
    
    # Impact parameter 
    b0 = (Gr*Per*Per*86400.*86400./(4*np.pi*np.pi))**(1/3.) * (MsMpR**(1/3.)) \
         * [(np.sin(phase*2*np.pi))**2. + (np.cos(inclin)*np.cos(phase*2*np.pi))**(2.)]**(0.5)

    # Model fit to data = light curve model * baseline flux * systematic model
    if transit==True:
        transit_curve=occultnl.occultnonlin(b0,np.sqrt(p[0]),p[18:22]) 
        model = transit_curve * p[1] * systematic_model_mpfit
    else:
        transit_curve=occultnl.occultnonlin(b0, rprs,p[18:22]) - 1.0
        eclipse=1.0 + p[0]*(1.0 + transit_curve/(np.max(transit_curve)-np.min(transit_curve)))
        model = eclipse * p[1] * systematic_model_mpfit

    resids = (y-model)/p[1]
    status=0
    
    return [status,(y-model)/err]

def systematic_model_grid_selection(size=4, transit=False):
    """ Returns model grid that indicates which parameters
    will be openly fit for and which will be fixed. Larger
    size will test higher powers. """
    
    if size not in [2,4,6]:
        sys.exit('Grid size can only be 2, 4, or 6')
    # MODEL GRID FOR JUST ECLIPSE DEPTH
    if size == 2:
        grid = [[0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], $ 
                [0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1], $ 
                [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1], $ 
                [0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1], $
                [0,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1]]
    if size == 4:
        # MODEL GRID UPTO THE 4th ORDER for HST & DELTA_lambda
        grid = [[0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

    if grid == 6:
        # MODEL GRID UPTO THE 6th ORDER for HST & DELTA_lambda
        grid = [[0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] ,$
                [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]

    if transit == True:
        grid[:,18:22]=0
        
    return grid

  
def whitelight_eclipse2017(p_start, img_date, allspec, plotting=False
                           , fixtime=False, norandomt=False, openinc=False, savefile=False
                           , transit=False):
    """
  NAME:                        
       WHITELIGHT_ECLIPSE2017.py
  
  AUTHOR:
     Based on Hannah R. Wakeford, NASA/GSFC code 693, Greenbelt, MD 20771
     hannah.wakeford@nasa.gov

     Kyle Sheppard: Converted to python and added eclipse functionality

  PURPOSE:
     Perform Levenberg-Marquardt least-squares minimization with
     MPFIT on spectral data from HST WFC3 to compute the band-integrated light curve

  MAJOR TOPICS:
     Generate band-integrated light curve
     Determine the most favoured systematic model for the observation
     Measure the marginalised secondary eclipse depth

  CALLING SEQUENCE:;
     whitelight_eclipse(p_start, img_date, allspec,plotting=False
                        fixtime=False, norandomt=False, openinc=False, savefile=False
                        , transit=False)

INPUTS:

 P_START - priors for each parameter used in the fit passed in
 an array in the form
 p_start = [rprs,epoch,inclin,MsMpR,Per,exposuretime, fp]

   fp - secondary eclipse depth  
   epoch - center of eclipse time
   inclin - inclination of the planetary orbit
   MsMpR - density of the system where MsMpR =
           (Ms+Mp)/(R*^3D0) this can also be calculated
           from the a/R* following
           constant1 = (G*Per*Per/(4*!pi*!pi))^(1D0/3D0) 
           MsMpR = (a_Rs/constant1)^3D0
           Per - Period of the planet in days
           exposuretime - exposure time for the total science exposure
   rprs - Planetary radius/stellar radius
   
 IMG_DATE - array of time of each exposure as extracted from the .fits header (MJD)


 ALLSPEC - 2D array each row containing the target stellar
           spectrum extracted from the exposure images.
           NOTE: this should be in units of e-/s when
           extracted from final science frame. If
           extracted by difference imaging this will be in
           units of e- see exposuretime settings for
           implications of this.

 SAVEFILE - root directory where you want to save the outputs
            e.g. '/Volumes/DATA1/user/HST/Planet/sav_file/'

 PLOTTING - set as True or False to see the plots

 FIXTIME - True to keep center of eclipse/transit time fixed
     
 NORANDOMT - True to not allow small random changes to center time

 OPENINC - True to fit for inclination

 TRANSIT - True for transit light curves, default false for eclipse 
 
    """
    
    # SET THE CONSTANTS 
    # constant = [GAIN,READNOISE,G,JD,DAY_TO_SEC,Rjup,Rsun,MJup,Msun,HST_SECOND,HST_PERIOD]
    constant = np.asarray([2.5,20.2,6.67259e-11,2400000.5,86400,7.15e7
                           ,6.96e8,1.9e27,1.99e30,5781.6,0.06691666])

    # TOTAL NUMBER OF EXPOSURES IN THE OBSERVATION
    nexposure = len(img_date)
    
    # CALCULATE THE SHIFT IN DELTA_lambda
    sh = np.zeros(nexposure)
    nLag=3
    for i in trange(nexposure, desc='Performing cross correlation'):
        # Subtract mean to mimic 
        inp1=allspec[-1,:]-allspec[-1,:].mean()
        inp2=allspec[i,:]-allspec[i,:].mean()
        corr_tuple = plt.xcorr(inp1, inp2,  maxlags=nLag)
        lag,corr=corr_tuple[0],corr_tuple[1]
        mx=np.argmax(corr)
        srad=3
        sublag=lag[max(mx-srad,0):max(mx+srad,2*nLag+1)]
        subcorr=corr[max(mx-srad,0):max(mx+srad,2*nLag+1)]
        p=np.polyfit(sublag, subcorr, 2)
        sh[i]=-p[1]/2./p[0] 


    # SET THE CONSTANTS USING THE PRIORS
    rprs = p_start[0] 
    epoch = p_start[1]
    inclin = p_start[2]
    MsMpR = p_start[3]
    Per = p_start[4]
    exposuretime = p_start[5]
    fp=p_start[6] #eclipse depth (planetary flux)
    flux0 = allspec[0,:].sum()
    T0 = img_date[0]

    if transit == False:
        depth=fp
    else:
        depth=rprs*rprs
    
    m = 0.0         # Linear Slope
    xshift1 = 0.0   # X-shift in wavelength
    xshift2 = 0.0   # X-shift^2 in wavelength
    xshift3 = 0.0   # X-shift^3 in wavelength
    xshift4 = 0.0   # X-shift^4 in wavelength
    xshift5 = 0.0   # X-shift^5 in wavelength
    xshift6 = 0.0   # X-shift^6 in wavelength
    HSTP1 = 0.0     # HST orbital phase
    HSTP2 = 0.0     # HST orbital phase^2
    HSTP3 = 0.0     # HST orbital phase^3
    HSTP4 = 0.0     # HST orbital phase^4
    HSTP5 = 0.0     # HST orbital phase^5
    HSTP6 = 0.0     # HST orbital phase^5
    c1 = 0.0        # Limb-darkening parameters
    c2 = 0.0 
    c3 = 0.0 
    c4 = 0.0 

    #PLACE ALL THE PRIORS IN AN ARRAY
    p0 = [depth,flux0,epoch,m,HSTP1,HSTP2,HSTP3,HSTP4,HSTP5,HSTP6,xshift1
          ,xshift2 ,xshift3,xshift4,xshift5,xshift6,inclin,MsMpR,c1,c2,c3,c4,Per,T0]
      
    nParam=len(p0)
    # SELECT THE SYSTEMATIC GRID OF MODELS TO USE 

    grid = np.transpose(systematic_model_grid_selection('four'), transit)
    # do I need to transpose?
    nsys = len(grid[:,0])

    #  SET UP THE ARRAYS  ;
    sys_depth = np.zeros(nsys,2)
    sys_model_x = np.zeros(nsys,4000)
    sys_model = np.zeros(nsys,4000) 
    sys_lightcurve_x = np.zeros(nsys,nexposure)
    sys_lightcurve = np.zeros(nsys,nexposure)
    sys_lightcurve_err = np.zeros(nsys,nexposure)
    sys_residuals = np.zeros(nsys,nexposure)
    sys_params = np.zeros(nsys,nParam)
    sys_params_err = np.zeros(nsys,nParam)
    sys_evidence = np.zeros(nsys) 
    sys_model_full=np.zeros(nsys,nexposure)

    phase = np.zeros(nexposure)
    HSTphase = np.zeros(nexposure)

    # Scatter of the residuals for each model
    resid_stddev = np.zeros(nsys)
  
    #run 1 AIC and parameters from the fit
    run1_AIC = np.zeros(nsys)
    run1_params = np.zeros(nsys,nParams)

    #  ITERATION RUN
    #  First run 4 trials with slightly shifted center of eclipse times
    #  and secondary eclipse depths 

    ntrials=5
    tcenter = np.zeros(ntrials+1)

    t1 = 5./60./24.
    tcenter[0] = epoch
    # HERE #
    if fixtime==False and norandomt==False:
        tcenter[1,:] = epoch + t1*np.random.normal(size=ntrials)

  
    # Test arrays ;
    AIC_test = np.zeros(ntrials+1)
    depth_test = np.zeros(ntrials+1)


    print '----------      ------------     ------------'
    print '          1ST FIT         '
    print '----------      ------------     ------------'

    x = img_date
    y=allspec.sum(axis=1) * exposuretime
    err = np.sqrt(y)
    phot_err=1e6/np.median(err)
    
    # Normalised Data
    # get in eclipse orbit, or first transit orbit
    ### Check if this works
    orbit_start, orbit_end=orbits('holder', x='x', y='y', transit=transit)[1]
    norm=np.median(y[orbit_start:orbit_end])
  
    rawerr=err
    rawflux=y
    err = err/norm
    y = y/norm
    flux0=1.0
    # Inflate noise more than photon uncertainty
    # err=err * 5.0   
 
    HSTphase = (img_date-T0)/constant(10)           
    HSTphase -= np.floor(HSTphase)
    HSTphase[HSTphase > 0.5] = HSTphase[HSTphase > 0.5] -1.0
  
    for s, systematics in tqdm(enumerate(grid), desc='First MPFIT run'):
        system=np.transpose(systematics)
        if fixtime==False and norandomt==False:
            for n in range(ntrials):
                print '-----------------'
                print 'TRIAL', n, '  Model:', s
                print 'Depth =', depth, '  Epoch = ',tcenter[n]
                print '  '

                epoch = tcenter[n]

                # Reset priors
                p0 = [depth,flux0,epoch,m,HSTP1,HSTP2,HSTP3,HSTP4,HSTP5,HSTP6,xshift1
                      ,xshift2,xshift3,xshift4,xshift5,xshift6,inclin,MsMpR,c1,c2,c3,c4,Per,T0]
                if openinc==True: systematics[16] = 0
                fa = {'x':x, 'y':y, 'err':err, 'sh':sh, 'rprs':rprs, 'transit':transit}
                parinfo=[]
                for i in range(len(p0)):
                    # Are value and system right?
                    parinfo.append({'value':p0[i], 'fixed':system[i], 'limited':[0,0],
                                    'limits':[0.,0.]})
    
                m = mpfit.mpfit(lightcurve,functkw=fa,parinfo=parinfo)
                params_test=m.params
                perror = m.perror
                nfree=len(systematics)-systematics.sum()
                dof=len(x)-nfree
                bestnorm=m.fnorm
                # For each tested epoch time, get the depth and the AIC
                # Make sure this is all computed correctly
                AIC_test[n] = bestnorm + nfree
                depth_test[n] = params_test[0]

        
            # Find the epoch time with the lowest AIC. Use it (or the average of 
            # values if multiple) to get starting depth and epoch time for
            # next round of fits.
            best = np.argmin(AIC_test)
            print 'Best eclipse depth prior =', fp_test[best]
            print 'Best center of eclipse prior =', tcenter[best]
            depth = np.mean(depth_test[best])
            epoch = np.median(tcenter[best])
  
    


        #Re-run the fitting process with the best prior as defined by the iterative fit
        p0 = [depth,flux0,epoch,m,HSTP1,HSTP2,HSTP3,HSTP4,HSTP5,HSTP6,xshift1
              ,xshift2,xshift3,xshift4,xshift5,xshift6,inclin,MsMpR,c1,c2,c3,c4,Per,T0]

        # MPFIT ;;;;;;;;;;;;;;;;;;;;;;;
        fa = {'x':x, 'y':y, 'err':err, 'sh':sh, 'rprs':rprs, 'transit':transit}
        if fixtime==True: systematics[2] = 1 
        if openi==True: systematics[16] = 0 
        parinfo=[]
        for i in range(len(p0)):
            parinfo.append({'value':p0[i], 'fixed':system[i], 'limited':[0,0],
                            'limits':[0.,0.]})
       
        m2 = mpfit.mpfit(lightcurve,functkw=fa,parinfo=parinfo)
        params_w=m2.params
        nfree=len(systematics)-systematics.sum()
        dof=len(x)-nfree
        perror=m.perror
        bestnorm=m.fnorm

        AIC = bestnorm + nfree
        pcerror = perror*np.sqrt(bestnorm/nfree)

        print 'Depth = ', params_w[0], ' at ', params_w[2]

        # Re-Calculate each of the arrays dependent on the output parameters
        phase = (x-params_w[2])/params_w[22] 
        phase -= np.floor(phase)
        phase[phase > 0.5] = phase[phase > 0.5] -1.0

        # ECLIPSE MODEL: calculate the eclipse model for the resolution of the data points
        # this routine is from MANDEL & AGOL (2002)

        # Calculate the impact parameter
        b0 = (constant[2]*params_w[22]*params_w[22]*86400.*86400./(4*np.pi*np.pi))**(1/3.) \
             * (params_w[17]*(1/3.)) * [(np.sin(phase*2*np.pi))**2.
                                        + (np.cos(params_w[16])*np.cos(phase*2*np.pi))**(2.)]**(0.5)

        systematic_model = ((phase*params_w[3] + 1.0) 
                            * (HSTphase*params_w[4] + HSTphase**2.*params_w[5]
                               + HSTphase**3.*params_w[6] + HSTphase**4.*params_w[7]
                               + HSTphase**5.*params_w[8] + HSTphase**6.*params_w[9] + 1.0)
                            * (sh*params_w[10] + sh**2.*params_w[11] + sh**3.*params_w[12]
                               + sh**4.*params_w[13]+ sh**5.*params_w[14] + sh**6.*params_w[15]
                               + 1.0))

        if transit==True:
            transit_curve=occultnl.occultnonlin(b0,np.sqrt(p[0]),p[18:22]) 
            w_model = transit_curve * p[1] * systematic_model_mpfit
        else:
            transit_curve=occultnl.occultnonlin(b0, rprs,p[18:22]) - 1.0
            eclipse=1.0 + p[0]*(1.0 + transit_curve/(np.max(transit_curve)-np.min(transit_curve)))
            w_model = eclipse * p[1] * systematic_model_mpfit

        w_corrected = y / (params[1] * systematic_model)   
        w_residuals = (y - w_model)/params[1]
        w_err = err/params[1]

        resid_stddev[s] = np.std(w_residuals)
        run1_AIC[s] = AIC
        run1_params[s,:] = params_w

    #######################################

    #Determine which of the systematic models initially gives the best fit
    top = np.argmin(run1_AIC)

    # Scale error by resid_stddev[top]
    std=resid_stddev[top]
    if np.median(err) < std:
        error=err*std/np.median(err)

    print '----------      ------------     ------------'
    print '         FINAL FIT        '
    print '----------      ------------     ------------'
    for s, systematics in tqdm(enumerate(grid), desc='Finat MPFIT run'):

        # Define the new priors as the parameters from the best fitting
        # systematic model

        p0=np.transpose(run1_params[s,:]) # transpose needed?

        fa = {'x':x, 'y':y, 'err':error, 'sh':sh, 'rprs':rprs, 'transit':transit}
        if fixtime==True: systematics[2] = 1 
        if openi==True: systematics[16] = 0 
        parinfo=[]
        for i in range(len(p0)):
            parinfo.append({'value':p0[i], 'fixed':system[i], 'limited':[0,0],
                            'limits':[0.,0.]})
       
        m2 = mpfit.mpfit(lightcurve,functkw=fa,parinfo=parinfo)
        params=m2.params
        nfree=len(systematics)-systematics.sum()
        dof=len(x)-nfree
        perror=m.perror
        bestnorm=m.fnorm
        AIC = bestnorm + nfree
        pcerror = perror*np.sqrt(bestnorm/nfree)

        print 'Depth = ', params[0], ' at ', params[2]

        # Re-Calculate each of the arrays dependent on the output parameters
        phase = (x-params[2])/params[22] 
        phase -= np.floor(phase)
        phase[phase > 0.5] = phase[phase > 0.5] -1.0

     
        # --------------------- #
        #        EVIDENCE       #
        sigma_points = np.median(error)
        Mpoint = systematics.sum()
        Npoint = len(x) 
        # EVIDENCE BASED ON AIC ;
        sys_evidence[s] = -Npoint*np.log(sigma_points)-0.5*Npoint*np.log(2*np.pi)-0.5*(AIC+nfree)
        
        b0 = (constant[2]*params[22]*params[22]*86400.*86400./(4*np.pi*np.pi))**(1/3.) \
             * (params[17]*(1/3.)) * [(np.sin(phase*2*np.pi))**2.
                                      + (np.cos(params[16])*np.cos(phase*2*np.pi))**(2.)]**(0.5)
        
        # copy? 

        systematic_model = ((phase*params_w[3] + 1.0) 
                            * (HSTphase*params_w[4] + HSTphase**2.*params_w[5]
                               + HSTphase**3.*params_w[6] + HSTphase**4.*params_w[7]
                               + HSTphase**5.*params_w[8] + HSTphase**6.*params_w[9] + 1.0)
                            * (sh*params_w[10] + sh**2.*params_w[11] + sh**3.*params_w[12]
                               + sh**4.*params_w[13]+ sh**5.*params_w[14] + sh**6.*params_w[15]
                               + 1.0))

        if transit==True:
            transit_curve=occultnl.occultnonlin(b0,np.sqrt(p[0]),p[18:22]) 
            w_model = transit_curve * p[1] * systematic_model_mpfit
        else:
            transit_curve=occultnl.occultnonlin(b0, rprs,p[18:22]) - 1.0
            eclipse=1.0 + p[0]*(1.0 + transit_curve/(np.max(transit_curve)-np.min(transit_curve)))
            w_model = eclipse * p[1] * systematic_model_mpfit

        w_corrected = y / (params[1] * systematic_model)   
        w_residuals = (y - w_model)/params[1]
        w_err = err/params[1]

        # Smooth Transit Model
        xsmooth = np.arange(4000)*0.00025-.5
        bsmooth = (constant[2]*params[22]*params[22]*86400.*86400./(4*np.pi*np.pi))**(1/3.) \
                  * (params[17]*(1/3.)) * [(np.sin(xsmooth*2*np.pi))**2.
                                           + (np.cos(params[16])*np.cos(xsmooth*2*np.pi))**(2.)]**(0.5)
        if transit==True:
            transit_curve=occultnl.occultnonlin(bsmooth,np.sqrt(p[0]),p[18:22]) 
            model_smooth = transit_curve * p[1] * systematic_model_mpfit
        else:
            transit_curve=occultnl.occultnonlin(bsmooth, rprs,p[18:22]) - 1.0
            eclipse=1.0 + p[0]*(1.0 + transit_curve/(np.max(transit_curve)-np.min(transit_curve)))
            model_smooth = eclipse * p[1] * systematic_model_mpfit
        
        # PLOTTING
        if ploton == True:
            plt.close()
            plt.errorbar(img_date, y, error,ecolor='red', color='red', marker='o')
            plt.ylim([0.982, 1.005])
            plt.plot(img_date, systematic_model, color='blue', marker='o')
            plt.errorbar(img_date, w_corrected, w_err, marker='x', color='green', ecolor='green')
            plt.show(block=False)

        # SAVE out the arrays for each systematic model ;
        sys_depth[s,0] = params[0]
        sys_depth[s,1] = pcerror[0]
        sys_lightcurve_x[s,:] = phase
        sys_lightcurve[s,:] = w_corrected
        sys_lightcurve_err[s,:] = w_err
        sys_model_x[s,:] = xsmooth
        sys_model[s,:] = model_smooth
        sys_residuals[s,:] = w_residuals
        sys_params[s,:] = params
        sys_params_err[s,:] = pcerror
        sys_model_full[s,:] = w_model

    #;;;;;;;;;;;;;;;;;;;
    #;;;;;;;;
    #;
    #;MARGINALIZATION!!!
    #;
    #;;;;;;;;
    #;;;;;;;;;;;;;;;;;;;

    # ------------------------------- ;
    #            EVIDENCE             ;
    aics = sys_evidence[;,1] 
    depth_array = sys_depth[:,0]         
    depth_err_array = sys_depth[:,1]     
    epoch_array = sys_params[:,2]       
    epoch_err_array = sys_params_err[:,2] 
    inc_array= sys_params[:,16]
    inc_err_array=sys_params_err[:,16]
    
    # Reverse sort as trying to MAXIMISE the negative log evidence
    a=np.sort(aics)[::-1] 
    best=np.argmax(aics)
    print best
    # print aics

    zero = np.where(aics < -500)
    if (len(zero) > 1): print 'Some bad fits - evidence becomes negative'
    if (len(zero) > 24):
        sys.exit('Over half the systematic models have nevative evidence, adjust and rerun')

    aics[aics < -500] = np.min(aics[aics>-500])

    beta=100.
    #beta = np.min(aics)
  
    w_q = (np.exp(aics-beta))/np.sum(np.exp(aics-beta))
    bestfit=np.argmax(w_q)
    n01 = np.where(w_q >= 0.1)
    
    stdResid = np.std(sys_residuals[bestfit,:]) 
    return_array[0]=stdResid
    print 'Evidences: ', aics
    print 'Weights: ', w_q

    print str(len(n01)) + ' models have a weight over 10%. Models: ', n01 , w_q(n01)
    print 'Most likely model is number ' +str(bestfit) +' at weight = ' + str(np.max(w_q))

    depth = depth_array
    depth_err = depth_err_array

    # Marganilze depth formula 15 and 16 from Wakeford 2016
    mean_depth=np.sum(w_q*depth)
    theta_qi=depth
    variance_theta_qi=depth_err*depth_err
    error_theta_i = np.sqrt(np.sum(w_q*((theta_qi - mean_depth)**2 + variance_theta_qi )))
    print 'Depth = %f  +/-  %f' % (mean_depth, error_theta_i)
    marg_depth = mean_depth
    marg_depth_err = error_theta_i 
    
    # Marganilze tcenter
    t0 = epoch_array
    t0_err = epoch_err_array
    inc=inc_array
    inc_err=inc_err_array
 
    print 'Depth'
    print depth[a]
    print depth_err_array[a]
    print 'Center Time'
    print t0[a]
    print t0_err[a]
    mean_t0=np.sum(w_q*t0)
    theta_t0=t0
    variance_theta_t0q = t0_err*t0_err
    error_theta_t0 = np.sqrt(np.sum(w_q*((theta_t0 - mean_t0)**2 + variance_theta_t0q )))
    print 'Central time = %f +/- %f' % (mean_t0, error_theta_t0) 
    marg_epoch = mean_t0
    marg_epoch_err = error_theta_t0
    
    # Marginalize over inclination
    if openi==True:
        mean_inc=np.sum(w_q*inc)
        bestfit_theta_inc=inc
        variance_theta_incq = inc_err*inc_err
        error_theta_inc = np.sqrt(np.sum(w_q*((bestfit_theta_inc - mean_inc )**2
                                              + variance_theta_incq )))
        print 'Inclination = %f +/- %f' % (mean_inc*180./np.pi, error_theta_inc)
        marg_inc_err = error_theta_inc
    else:
        marg_inc=inc[0]
        marg_inc_err=0
 
    plt.errorbar(sys_lightcurve_x[bestfit,:], sys_lightcurve[bestfit,:], sys_lightcurve_err[bestfit,:]
                 ,marker='o', color='b', ecolor='b')
    plt.plot(sys_model_x[bestfit,:], sys_model[bestfit,:], marker='^')
    plt.show()

    rms=np.std(sys_residuals[bestfit,:])*1e6
    ratio=rms/phot_err
    print rms, phot_err, ratio
    print marg_depth, marg_depth_err

    if savewl:
        #################3 make sure this works
        ### Two dataframes, both multi-indexed
        # First: Anything dependent on systematic models. Adjusted data, weights,
        # residuals, etc.
        # Second: All data things for plots. Flux/errors normalized and un. Dates,
        # And then all relevant fitting results (RMS/photon)

        # To retrieve: df.loc[visit, type you want][column]
        # Example: wl_models_info.loc['hatp41/visit01','Params']['Model 12'].values[0]
        
        cols = ['Model ' + str(i) for i in range(nsys)]
        subindex=['Weight'] + ['Corrected Flux']*nexposure + ['Corrected Phase']*nexposure
                  + ['Corrected Error']*nexposure + ['Residuals']*nexposure
                  + ['Params']*nParam + ['Params Errors']*nParam + ['AIC Evidence']
                  + ['Smooth Model']*4000 + ['Smooth Model Phase']*4000
        ind=pd.MultiIndex.from_product([[savewl], subindex])
        wl=pd.DataFrame(np.vstack(w_q, sys_lightcurve, sys_lightcurve_x
                                  , sys_lightcurve_err, sys_residuals
                                  , sys_params, sys_params_err, sys_evidence
                                  , sys_model,sys_model_x), columns=cols, index=ind)

        ind2a=pd.MultIndex.from_product([[savewl],['data']*nexposure])
        colsa=['Obs Date', 'Normalized Flux', 'Flux', 'Normalized Error'
              , 'Error', 'Wavelength Shift']
        dataa=np.vstack(img_date, y, rawflux, error, rawerr, sh)
        colsb='Values'
        datab=[marg_depth, marg_depth_err, marg_epoch, marg_epoch_err, rms
               , phot_err, ratio, orbit_start, orbit_end, rprs, stdResid]
        ind2b=pd.MultiIndex.from_product([[savewl],['Marg Depth', 'Depth err'
                                                  , 'Marg Epoch', 'Epoch err', 'RMS'
                                                    , 'photon err' , 'ratio', 'Norm index1'
                                                    , 'Norm index2', 'RpRs', 'Std Best Resids']])
        df1 = pd.DataFrame(dataa, columns=colsa, index=ind2a)
        df2 = pd.DataFrame(datab, columns=colsb, index=ind2b)
        wl_data = pd.concat((df1,df2))
    
        try:
            cur=pd.read_csv('./wl_models_info.csv')
            cur=pd.concat((cur,wl_smooth))
            cur.to_csv('./wl_models_info.csv', index=False)
        except IOError:
            wl.to_csv('./wl_models_info.csv', index=False)
            
        try:
            curr=pd.read_csv('./wl_data.csv')
            curr=pd.concat((curr,wl_data))
            curr.to_csv('./wl_data.csv', index=False)
        except IOError:
            wl_data.to_csv('./wl_data.csv', index=False)

    return np.asarray([marg_depth, marg_depth_err, marg_epoch
                       , marg_epoch_err, marg_inc, marg_inc_err])




