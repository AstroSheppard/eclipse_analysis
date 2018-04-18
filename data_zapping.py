import glob
import sys
import os
import shutil

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm 

# Basically debugged. Make sure median pixel zapping works for actual data

coords=[]

def onclick(event):
    global ix, iy
    ix=event.xdata
    iy=event.ydata
    print 'x = %d, y = %d'%(
        ix, iy)

    global coords
    coords.append((ix, iy))
    
    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)

        
def test(img, window, final_img):
        """ Test to make sure window size does not cut 
        out important data"""
        scale=46.7
        raw_fit=bkg(img, window, test=True)
        final_img=final_img*scale
	fits.writeto('./test_zapping/raw.fits', final_img, overwrite=True)
	sys.exit('Compare background-subtracted image with input image (raw.fits)')

def get_data(visit):
        """ Extract only quality, spectroscopic fits files """
        # Read in data and decare arrays
        raw=np.asarray(glob.glob('./'+visit+'/*ima.fits'))
        fit=fits.open(raw[0])
        example=fit[1].data
        xsize, ysize=np.shape(example) # make sure x and y are correct length
        fit.close()

	# Make sure the data is spectroscopic and SPARS10 (scan, not rapid)
        obstype=np.zeros(len(raw)).astype(str)
        rate=np.zeros(len(raw)).astype(str)
        quality=np.zeros(len(raw)).astype(str)
  	for i,img in tqdm(enumerate(raw)):
                fit=fits.open(img)
                header=fit[0].header
      	        obstype[i]=header['OBSTYPE'] 
                rate[i]=header['SAMP_SEQ']
            	quality[i]=header['QUALITY']
                fit.close()
        index=(obstype == 'SPECTROSCOPIC')* (
            (rate ==  'SPARS10') + (rate == 'SPARS25')) * (quality != 'LOCKLOST')
        og=len(raw)
        raw=raw[index]
        print len(raw), " images remain out of", og, "originals"
        return raw
    
def zapped(allspec):
        """Input is 3D numpy array of all bkg-removed exposure. 
        Median values and sigma cuts are used to remove cosmic rays"""
        dims=np.empty_like(allspec)
	allspecout=dims
        aspec=dims
  	med_rows=dims
        nspec, nx, ny = allspec.shape
  	nloop = 2                      
  	nzap = 0                          
  	for j in range(nloop):      
     		if j==0:
                        n=8
                else:
                        n=6

                # Normalize each row by it's median to account
                # for uneven scan rates. Loop over pixels, and
                # check for cosmic rays by comparing each pixel to
                # its median value over all exposures and finding
                # stddev. If stddev is above a threshold, then set that
                # to median of it's time series value. Make sure
                # to scale back to pre-row normalization
                
                rows=np.median(allspec, axis=1)          
                np.place(rows,rows==0,1)
                #rows[z,y] gives median value of row at column y in image z
		for x in range(nx):
                        aspec[:,x,:]=allspec[:,x,:]/rows
			for y in range(ny): 
				# Find median of every pixel's time series
				med = np.median(aspec[:,x, y]) 
				sigma = np.sqrt(np.sum((aspec[:,x, y]-
                                                        med)**2)/nspec)   
				for z in range(nspec): 
					dif = np.abs(aspec[z,x,y]-med)
					if dif > n*sigma:
                                                allspec[z,x,y] = med*rows[z,y]
        
     		print 'Number Zapped is ' + str(nzap)
  	return allspec

def maxrow(frame):
        """Return index of maximum flux row"""
	f=np.sum(frame,axis=1)
  	f1=np.argmax(f)	
  	return f1

# Input 1 fits file, make array with each extension (3-D: [ext, x, y])
def bkg(raw, size_window, test=False):
        """ Use differential frames (Deming 2013) to remove
        background."""
        with fits.open(raw) as exp:
                xlen, ylen = np.shape(exp[1].data)
                header=exp[0].header
                corr=header['UNITCORR']
                nFrames=exp[-1].header['EXTVER']
        
  	        frames = np.zeros((nFrames, xlen, ylen))
  	        times=np.zeros(nFrames)
  	        frame_diffs = np.zeros((nFrames-1, xlen, ylen))
  	        count=0
                # Iterate through extensions in fits file
  	        for item in exp[1:]:
                    if 'SCI' in item.header['EXTNAME']:
                        frames[count,:,:]=item.data
     		        times[count]=item.header['SAMPTIME']
     		        count= count + 1
  	

	# Check if the units are in electrons or in electrons/s. If per
	# second, multiply frame by samptime. 
	# Zero out background that is far enough from source to be noise (or
	# is from secondary source)
  	if 'OMIT' in corr:
     		# Last frame is for tsamp=0. If brightness not in units of e/s, then
		# this last frame will give large negatives. Set it to 0

                # Adjust window size here, since "dir" is not a perfect correlation.
     		window=np.floor(size_window*2.4).astype(int)
     		frames[-1,:,:]=0.0
     		ny=frames.shape[2]
    		# window=40
    		# window=12 ; hatp41
     		for j in range(count-1):
        		f1=frames[j,:,:]      
        		f2=frames[j+1,:,:]     
        		frame_diffs[j,:,:]=f1-f2
        		mrow=maxrow(frame_diffs[j,:,:])
                       
                        # prevent window outside of data
                        if mrow+window > ny-1: window = ny-1-mrow
                        # Mask data, find median of background,
                        # subtract from image, zero out bkg
                        bg=np.concatenate((frame_diffs[j,:,:mrow-window],
                                           frame_diffs[j,:,mrow+window:-1]),axis=1)
                        med=np.median(bg)

                      
       			frame_diffs[j,:,:]=frame_diffs[j,:,:]-med
        		frame_diffs[j,:mrow-window,:]=0 
        		frame_diffs[j,mrow+window:,:]=0
                        if test:
        	                file1='./test_zapping/dif'+str(j)+'.fits'
        	                file2='./test_zapping/frame'+str(j)+'.fits'
                                fits.writeto(file1,frame_diffs[j,:,:], overwrite=True)
                                fits.writeto(file2,f1, overwrite=True)

  	if 'COMPLETE' in corr:
    		window=size_window
     	        ny=frames.shape[2]
     		for j in range(count-1):
        		f1=frames[j,:,:]*times[j]
        		f2=frames[j+1,:,:]*times[j+1]
        		frame_diffs[j,:,:]=f1-f2
        		mrow=maxrow(frame_diffs[j,:,:]) 
        		if mrow+window > ny-1: window = ny-1-mrow 
        		# Find median of background noise
                        bg=np.concatenate((frame_diffs[j,:,:mrow-window],
                                           frame_diffs[j,:,mrow+window:-1]))
                        med=np.median(bg)
       			frame_diffs[j,:,:]=frame_diffs[j,:,:]-med
        		frame_diffs[j,:mrow-window,:]=0 
        		frame_diffs[j,mrow+window:,:]=0
                        
                        if test:
        	                file1='./test_zapping/dif'+str(j)+'.fits'
        	                file2='./test_zapping/frame'+str(j)+'.fits'
                                fits.writeto(file1,frame_diffs[j,:,:], overwrite=True)
                                fits.writeto(file2,f1, overwrite=True)
  
 
  
  	output=frame_diffs.sum(0)
        if test: fits.writeto('./test_zapping/fbkg.fits', output, overwrite=True)
        return output
  

if len(sys.argv) < 3:
    sys.exit('Please use python [data_zapping.py] [planet] [visit]')
visit=sys.argv[1]+'/'+sys.argv[2]
raw=get_data(visit)
n_forward, n_reverse=0,0
direction=np.zeros(len(raw))
print direction
for i,img in tqdm(enumerate(raw), desc='Getting data'):
    exp=fits.open(img)
    header=exp[0].header
    dire=header['POSTARG2']
    direction[i]=dire
    if dire > 0:
       	n_forward+=1
        if n_forward == 1:
            forward_img=exp[1].data
            exp.close()
    else:
        n_reverse+=1
       	if n_reverse == 1: 
            reverse_img=exp[1].data
            exp.close()
# Center the source                               
if n_forward == 0:
    img=reverse_img
else:
    img=forward_img            
fig=plt.figure()
ax=plt.imshow(img)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
print "Click the top-left then the bottom-right corners"
plt.show()
# need to make sure x and y are being extracted correctly
coords= [int(i) for item in coords for i in item]
x1=coords[1]
x2=coords[3]
y1=coords[0]
y2=coords[2]
xlen=x2-x1
ylen=y2-y1

if n_forward > 0:
    allimages_f = np.zeros((n_forward, xlen, ylen))
    allheader_f = np.asarray([])
if n_reverse > 0:
    allimages_r = np.zeros((n_reverse, xlen, ylen))
    allheader_r = np.asarray([])

# determine window size for both forward and reverse images
# direction array contains something related to scan rate for each image
w1=direction[0]
w2=direction[1]
rwindow=0
fwindow=0
        
# If direction[0] is negative, then it is a reverse scan, and
# we set the window to be a size that typically captures all
# source photons
if w1 < 0:
    rwindow=np.ceil(2*np.abs(w1))
# If it's positive, then it's a forward scan and we set the
# forward window instead
else:
    fwindow=np.ceil(2*np.abs(w1))
# Now we check the second scan, which can be negative for a
# bi-directional, in which case we set the reverse window. 
if w2 < 0:
    rwindow=np.ceil(2*np.abs(w2))
# For uni-directional, the rate doesnt change so this does nothing
else:
    fwindow=np.ceil(2*np.abs(w2))
# For one data set the scan direction (I use as a proxy for rate)
# was super low, so we can set it to be either the other directions
# window or 1, in case the other direction is 0 (ie,
# unidirectional and low scan "rate")
if rwindow <= 1: rwindow = max(fwindow,1)
if fwindow <= 1: fwindow = max(rwindow,1)

if len(sys.argv)==4: test(raw[1], rwindow, reverse_img) # ensure window size is okay

# For each exposure, check if reverse or forward scan.
# Subtract background
# Apply window to center source
# Save image to all images, save header to all headers 
f=0
r=0
for i,expo in tqdm(enumerate(raw), desc='Subtracting background'):
    hdr=fits.open(expo)[0].header.tostring(sep='\\n')
    if direction[i] > 0:
        raw_fit=bkg(expo, fwindow)
        # make sure x,y indexing is right. Same with appending/header
	img = raw_fit[x1:x2,y1:y2]
	allimages_f[f,:,:]=img
	allheader_f=np.append(allheader_f,hdr)
	f+=1
    else:
        raw_fit=bkg(expo, rwindow)
	img = raw_fit[x1:x2,y1:y2]
	allimages_r[r,:,:]=img
        allheader_r=np.append(allheader_r,hdr)
	r+=1
                
# Make directories for cleaned data/clear current directories

if n_reverse > 0:
    zapped_r=zapped(allimages_r)
    rzdir='./zapped2017/'+visit+'/reverse/cr/'
    r2dir='./zapped2017/'+visit+'/reverse/final/'
    try:
        os.makedirs(rzdir)
    except OSError:
        if os.path.isdir(rzdir):
            shutil.rmtree(rzdir)
            os.makedirs(rzdir)
        else:
            raise
    try:
        os.makedirs(r2dir)
    except OSError:
        if os.path.isdir(r2dir):
            shutil.rmtree(r2dir)
            os.makedirs(r2dir)
        else:
            raise
if n_forward > 0:
    zapped_f=zapped(allimages_f)
    fzdir='./zapped2017/'+visit+'/forward/cr/'
    f2dir='./zapped2017/'+visit+'/forward/final/'
    try:
        os.mkdir(fzdir)
    except OSError:
        if os.path.isdir(fzdir):
            shutil.rmtree(fzdir)
            os.mkdir(fzdir)
        else:
            raise
    try:
        os.mkdir(r2dir)
    except OSError:
        if os.path.isdir(f2dir):
            shutil.rmtree(f2dir)
            os.mkdir(r2dir)
        else:
            raise

# Write reduced data to a directory
for k in range(n_forward):
    filename = fzdir + "%03d"%k + 'f.zap.fits'
    zapp_image = zapped_f[k,:,:]
    hdr=fits.Header.fromstring(allheader_f[k], sep='\\n')
    fits.writeto(filename, zapp_image, header=hdr, overwrite=True)
for k in range(n_reverse):
    filename = rzdir + "%03d"%k + 'r.zap.fits'
    zapp_image = zapped_r[k,:,:]
    hdr=fits.Header.fromstring(allheader_r[k], sep='\\n')
    fits.writeto(filename, zapp_image, header=hdr)    
print 'Finished zapping ' + visit + ' visit'


  
 


