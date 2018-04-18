import glob
import sys

import numpy as np
from astropy.io import fits
from tqdm import tqdm
from tqdm import trange
from scipy.special import erfinv
import matplotlib.pyplot as plt


def bad_pixels(visit, direction):
    folder='./zapped2017/'+visit+'/'+direction+'/cr/*.fits'
    data=glob.glob(folder)
    n=len(data)
    x=np.zeros(3)                   #[apt,min,max]
    y=np.zeros(3) 
    xl=np.zeros(3) 
    yl=np.zeros(3) 
    if n > 1:
        ####Find base aperture####
        x,y=aperture(data[1])
        # Use larger aperture to include all light
        update=np.asarray([1,-1,1])
        xl=x+4*update
        yl=y+4*update
        for i,img in tqdm(enumerate(data), desc='Correcting bad pixels'):
            exposure=fits.open(img)
            print len(exposure)
            img_data=exposure[0].data #### is this 0 right?
            expo=img_data[xl[1]:xl[2], yl[1]:yl[2]]
            hdr=exposure[0].header
            exposure.close()
            ####
            zapped=pixel_zapping(xl[0], yl[0], expo)
            # save zapped images in folder named by dataset
            filename = './zapped2017/'+ visit + '/'+direction+'/final/' + "%03d"%i +'.zap.fits'
            fits.writeto(filename,zapped, header=hdr, overwrite=True)

def aperture(data):
    """Given a sample exposure, return the basic aperture."""
    # Read in list of fits files, return the x and y apertures
    exp = fits.open(data)
    exposure=exp[0].data
    xlen,ylen = exposure.shape
    exp.close()
    center=np.zeros(2).astype(int)
    exposure_xcen = exposure[(xlen/2),:]
    
    # Use 10% of max pixel as gauge for edge cutoff
    max_count=np.max(exposure)
    scan = np.where(exposure_xcen > max_count/10.)[0]
    scan_width = len(scan) 
    scan_center = int(np.median(scan))

    exposure_ycen = exposure[:,scan_center] 
    height = np.where(exposure_ycen > max_count/10.)[0]
    scan_height = len(height)
    height_cen=int(np.median(height))

    # Redo width for new y-center to make sure
    # center is accurate
    exposure_xcen=exposure[height_cen,:]
    scan2=np.where(exposure_xcen > max_count/10.)[0]
    scan_width= len(scan2)
    scan_center=int(np.median(scan2))
    center[0]=height_cen
    center[1]=scan_center
    xapt=scan_height/2
    yapt=scan_width/2
    xmin=center[0]-xapt
    xmax=center[0]+xapt
    ymin=center[1]-yapt
    ymax=center[1]+yapt
    x=np.asarray([xapt,xmin,xmax])
    y=np.asarray([yapt,ymin,ymax])
    return (np.vstack((x,y)))


def pixel_zapping(xapt, yapt, allspec, plot=False):
    """ Removes bad pixels by comparing each pixel
    to its column's median."""
    nloop=3                       
    # Pixels expected to be outside sigma factor in image
    nextra=5.   
    np.place(allspec, allspec<0, 0)
    print 'Negative pixels = ', np.sum(allspec<0)

    # Zoom in to ignore lower-flux pixels in median calculation
    xapt1=xapt-4 
   

    inp=1-nextra/4./yapt/xapt1
    n=(2.)**(.5)*erfinv(inp)
    # Sigma rejection factor chosen such that only 5 pixels in image can be expected
    # to naturally lie outside the sigma range. I then can correct the pixels
    # outside the sigma range without worry about overcorrecting

    for j in trange(nloop):       
        # Take care of any other bad pixels by comparing each pixel to the
        # median of it's own column
        allspec1=allspec[4:-4,:]
        column_med=np.median(allspec1,axis=0)
        col_med=np.broadcast_to(column_med, (2*xapt, 2*yapt))
        column_med=np.broadcast_to(column_med, (2*xapt1, 2*yapt))
        sigma=np.sqrt(np.sum((allspec1-column_med)**2,axis=0)/(2*xapt1))
        sigma=np.broadcast_to(sigma, (2*xapt, 2*yapt))
        dif=np.abs(allspec-col_med)
        # Account for edge effects: If the median of the row is much much
        # different then the median of all rows, then don't change the pixels
        row_m=np.broadcast_to(np.median(allspec, axis=1), (2*yapt, 2*xapt)).T
        ycut=np.median(row_m)*.9
        yhigh=np.median(row_m)*1.1
        # Replace the bad pixels with the median of their column
        index=(dif > n*sigma)* (row_m >= ycut) * (row_m <= yhigh)
     
        allspec[index]=col_med[index]
        print 'Bad pixels found =',np.sum(index)

    zapped=allspec
    return zapped


if len(sys.argv) != 3:
    sys.exit('Use python2 [fullzap.py] [planet] [visit]')
else:
    visit=sys.argv[1]+'/'+sys.argv[2]
    bad_pixels(visit, 'forward')
    bad_pixels(visit, 'reverse')
