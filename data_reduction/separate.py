import sys


import numpy as np
import glob
from astropy.io import fits
from shutil import copyfile

def separate(planet):

    """ Go to all data file in a directory and separate them based on date
    (so that different visits are separated)"""

    # Read in files
    # Find the date of each
    # if date2-date1 greater than a hubble orbit period * 2, then save
    # files in new directory

    data=np.asarray(glob.glob('../planets/%s/*ima.fits' % planet))
    date=np.zeros_like(data)
 
    ; Read in dates for each exposure
    for i, img in enumerate(data):
        fit=fits.open(img)
        date[i]=(fit[0].header['EXPSTART']+fit[0].header['EXPEND'])/2.
        fit.close()

    time=np.zeros(len(data)-1)
    visit=np.zeros_like(data)
    hst_period=95.47
    for i in nspec-2:
        t=np.abs(date[i+1]-date[i])
        time[i]=t*24*60
        
        #  If time between data points is greater than 
        #  3 HST orbits from previous exposure,
        #  then I classify it as a new observation
        
        t=t/hst_period
        if t > 3: visit[i]=1

 
    nObs=np.sum(visit)+1
    fNames=np.arange(nObs)
   
    dirs=['../planets/%s/visit%02i/' % (planet, name) for name in fNames]
  
    for dir in dirs:
        try:
            os.makedirs(dir)
        except OSError:
            if os.path.isdir(dir):
                shutil.rmtree(dir)
                os.makedirs(dir)
            else:
                raise
    
    direct=dirs[0]
    numV=0
    for i, img in enumerate(data):
        copyfile('../%s/%s' % (planet, img), direct)
        if visit[i] == 1:
            numV=numV+1
            direct=dirs[numV]
   
 
