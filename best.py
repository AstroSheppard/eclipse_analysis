import os
import glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def wmean(val, dval):
   dw=np.sum(1.0/dval**2)   
   sum=np.sum(1.0*val/dval**2)/dw
   error=np.sqrt(1.0/dw)
   return [sum,error]

def best(planet):

   dir='./spectra_april/'
   data=glob.glob(dir+planet+'*.sav')
   # glob may be unnecessary if we save in one csv instead of many sav's
###   restore, data[0]
### probably will replce by saving variables as pandas DF in csvs instead of save files
   nBin=len(depth)
   nVisit=len(data)
   margdepths=np.zeros(nVisit)
   margerrs=np.zeros(nVisit)
   spectra=np.zeros((nVisit, nBin, 2))
   basedir='./WLresults/'
   for i in range(len(nVisit)):

  ###    restore, data[i]
      spectra[i,:,0]=depth  #? make sure indexing is correct
      spectra[i,:,1]=error
      if i == 0: wcenter=center
         WLfile=basedir+os.path.basename(data[i])
  ###    restore, WLfile
      margdepths[i]=marg_depth
      margerrs[i]=marg_depth_err

   
   margWL=wmean(margdepths, margerrs)
   dif=(-1.0*margdepths+margWL[0])*1e6

   for i in range(len(nVisit)):
      spectra[i,:,0]=spectra[i,:,0]+dif[i]

   final_spectra=np.zeros((nBin,2))
   for i in range(len(nBin)):
     ### bindepth=reform(spectra[:,i,0])
      bindepth=spectra[:,i,0].flatten() #? maybe flatten(1), maybe unnecessary. get shape
      binerror=spectra[:,i,1].flatten()
      final_bin=wmean(bindepth, binerror)
      final_spectra[i,:]=final_bin
   
   # Plot spectra
   center=wcenter
   depth, error=final_spectra[:,0], final_spectra[:,1]

   first=0
   last=-1
   center=center[first:last]  #? indexing the same? in python this cuts last element
   depth=depth[first:last]
   error=error[first:last]
   plt.errorbar(wcenter, depth, error, ls='', marker='o', color='red', ecolor='red')
   plt.xlabel('Wavelength')
   plt.ylabel('Eclipse Depth')
   plt.title('Emission Spectra')
   plt.savefig(dir+planet+'.pdf')
   plt.show()  
#   plt.close()

###   SAVE, filename=dir+planet+'.sav', center, results, depth, error, dif
# how to make this a pandas df? planet (repeated), center, depth, error, ?results, 


  




