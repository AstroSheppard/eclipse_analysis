import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import occultnl.occultnonlin as occultnl
import whitelight2018 as wl

def eclipse_time(date, properties):
    """Program to determine the expected eclipse time
     Inputs
     date: 1D array of the date of each exposure (MJD)
     properties: 1D array containing the last observed eclipse 
     and the period. (MJD, days)"""
    time=properties[1]
    period=properties[4]
    while time < date[0]:
        time+=period
    return float(time)

def get_orbits(date):
    """Procedure to organize light curve data by HST orbit"""
    orbit=np.zeros(1)
    for i=0 in range(len(date)-1):
        t=date[i+1]-date[i]
        if t*86400 > 1200.:
            orbit=np.append(orbit, i+1) # 1800s is about half an HST orbit
    return np.append(orbit, len(date))

def inputs(data, transit=True):

    """ Function to read in priors for a system. 
    INPUTS:
    data: data table of priors for a particular planet
    OUTPUTS:
    Returns array of system properties: [rprs^2, central event time, inc
    ,mpmsr, period, 1, eclipse depth] 
    """
    inp_values=pd.read_table(data,sep=' ')
    data_arr=inp_values.iloc[:,2].values
    labels=inp_values.iloc[:,0].values
    # Make sure format of table is right. Rn, label - units - values 
    ncols=3
  

    # Rj-m, Rsolar-m,deg-radian,days-secs,AU-m
    conversions=np.array([6.9911e7,6.957e8,np.pi/180.,86400,1.49598e11])
    G=6.67408e-11
    inc=data_arr[5]*conversions[2]
    period=data_arr[4]
    if transit==True and 'transit' in labels[6] or transit==False and 'transit' not in labels[6]:
        epoch=data_arr[6]-2400000.5
    else:
        epoch=data_arr[6]-2400000.5+period/2.
   
    a_R=data_arr[7]*conversions[4]/(data_arr[1]*conversions[1])
    MpMsR=(a_R**3)*4*np.pi*np.pi/(G*period*period*conversions[3]*conversions[3])
    rprs = data_arr[0]*conversions[0]/(data_arr[1]*conversions[1])
    depth = rprs*rprs*(data_arr[2]/data_arr[3])
  
    props=np.zeros(7)
    props[0]=rprs*rprs
    props[1]=epoch
    props[2]=inc
    props[3]=MpMsR
    props[4]=period
    props[5]=1
    props[6]=depth/3 # Eclipse depth prior: why div 3?
    return props

def correction(inputs, date, flux, transit=False):
    t0=inputs[1]
    depth=inputs[6]
    MsMpR=inputs[3]
    inclin=inputs[2]
    Per = inputs[4] 
    rprs=np.sqrt(inputs[0])
    JD = 2400000.5              
    Gr = 6.67259e-11
    
    phase = (date-t0)/Per
    phase = phase - np.floor(phase)
    phase[phase > 0.5] = phase[phase > 0.5] - 1.0
 
    
    b0 = (Gr*Per*Per*86400.*86400./(4*np.pi*np.pi))**(1./3) * (MsMpR**(1./3)) \
         * [(np.sin(phase*2*np.pi))**2. + (np.cos(inclin)*np.cos(phase*2*np.pi))**(2.)]**(0.5)

    if transit==True:
        model=occultnl(b0,rprs,[0,0,0,0]) 
    else:
        transit_curve=occultnl(b0, rprs,[0,0,0,0]) - 1.0
        model=1.0 + depth*(1.0 + transit_curve/(np.max(transit_curve)-np.min(transit_curve)))
        
    corrected=flux/model
    return corrected


def remove_bad_data(light_curve, spectra, light_corrected, date1, check=False, user_inputs):
    """Procedure to remove "bad" data from light curve"""
 
    med= np.median(light_corrected)
    sigma = np.sqrt(np.sum((light_corrected-med)**2)/(2*len(light_corrected)))
    medi=np.zeros_like(date1)+med

    sig3=medi+3*sigma
    sig4=medi+4*sigma
    sig5=medi+5*sigma
    sig3m=medi-3*sigma
    sig4m=medi-4*sigma
    sig5m=medi-5*sigma

    if check==False:
        nPasses=user_inputs[3]
        sigma_cut_factor=user_inputs[2]

    else:
        data=plt.plot(date1, light_corrected, ,'bo', xlabel='MJD', ylabel='Total Flux',ls='dotted')
        s5=plt.plot(date1, sig5,'pink',sig5m, 'pink')
        s5[0].set_label('5-sigma')
        s4=plt.plot(date1, sig4,'g', sig4m, 'g')
        s4[0].set_label('4-sigma')
        s3=plt.plot(date1, sig3,'r', sig3m, 'r')
        s3[0].set_label('3-sigma')
        plt.plot(date1, medi, label='Median',ls='solid')
        plt.legend(scatterpoints=1)
        plt.show(block=False)
        cut = raw_input("Enter the sigma-cut factor (3-5 recommended): ")
        sigma_cut_factor = float(cut)
        user_inputs[2]=sigma_cut_factor
        passes=raw_input("Enter the number of passes for the sigma-cut: ")
        nPasses=int(passes)
        user_inputs[3]=nPasses
        plt.close()
  
    # Cut out the "bad" data

    for j in range(nPasses):
        med= np.median(light_corrected)
        sigma = np.sqrt(np.sum((light_corrected-med)**2)/(2*len(light_corrected)))
        dif= np.abs(light_corrected-med)
        index=np.where(dif > sigma_cut_factor*sigma)
        light_curve=light_curve[index]
        date1=date1[index]
        light_corrected=light_corrected[index]
        spectra=spectra[index,:]
        ### Does this remove data correctly?
    return [light_curve, spectra, light_corrected, date1]


def preprocess_whitelight(visit, direction, x=0, y=0, ploton=True
                          , check=True, inp_file=False, savedata=False
                          , transit=False):

    """ 
    PURPOSE: Allow user to e xtract relevant orbital data from reduced time 
    series of a visit. Also allow user to exclude any outlier data points. 
    
    INPUTS

     x, y, allow the user to reduce aperture
     checks: set to "on" to manually reduce data

     If checks is set to on, "user_inputs" will return the inputs
     that the user used: [first orbit, last orbit, sigma cut factor, 
     number of passes, center eclipse time]. If checks is set to off, then
     the user_inputs array will be used as inputs (easier to automate) """
 
    folder = './zapped2017/%s/%s/final/*.zap.fits' % (visit, direction)
    data=np.asarray(glob.glob(folder))
    nexposure = len(data)
    print 'There are %d exposures in this visit' % nexposure

    date=np.zeros_like(data)
    time=np.zeros_like(data) # comment out
    
    test=fits.open(data[0])
    xlen, ylen = test[1].data.shape
    test.close()
    xlen-=2*x
    ylen-=2*y
    allspec=np.zeros((len(data),xlen, ylen))

    xmin=x
    xmax=xlen-x
    ymin=y
    ymax=ylen-y

    for i, img in enumerate(data):
        expfile=fits.open(img)
        hdr=exp[0].header
        exp=exp[1].data
        expfile.close() # this might need to be after below stuff
        date[i]=(hdr['EXPSTART']+hdr['EXPEND'])/2.
        time[i]=hdr['EXPTIME']
        expo=exp[xmin:xmax, ymin:ymax]
        allspec[i,:,:]=expo

    props=inputs('../planet/%s/inputs.dat' % visit, transit)
    orbit = np.zeros(1)
    folder='./' 

    date_order=np.argsort(date)
    date=date[date_order]
    allspec=allspec[date_order,:,:]

    # Classify the data by each HST orbit. Returns array (orbit) 
    # which contains the indeces for the start of each orbit

    orbit=get_orbits(date)

    print "Number of total orbits: %d" % len(orbit)-1

    # Choose which orbits to include in the eclipse fitting. 1-2 on either
    # side of the eclipse is recommended

    if check == False:
        if inp_file == True:
            df=pd.read_csv('./preprocess_info.csv')
            df=df[df.loc[:,'Transit']==transit]
            user_inputs=df.loc[visit+direction,'User Inputs'].values
        else: 
            sys.exit('Either allow checking or give csv file with pd info.')

        first_orbit=user_inputs[0]
        last_orbit=user_inputs[1]
        date1=date[orbit[first_orbit]:orbit[last_orbit+1]]
        allspec1=allspec[orbit[first_orbit]:orbit[last_orbit+1],:,:] 
        allspec1=np.sum(allspec1,axis=2) #spectra for each exposure: these axes may be backwards
        light = np.sum(allspec1, axis=1) # total light for each exposure

    if check == True:
        user_inputs=np.zeros(5)
        while check==True:
            if ploton==True:
                plt.plot(date, np.sum(allspec, (1,2)),'o')
                plt.xlabel('MJD')
                plt.ylabel('Total Flux')
                plt.show(block=False)
            first = raw_input("Enter the first orbit to include (starting from 0): ")
            first_orbit=int(first) # fix? also how to read prompts?
            user_inputs[0]=first_orbit
            last= raw_input("Enter the last orbit to include (starting form 0): ")
            last_orbit=int(last)
            if ploton==True: plt.close()
            user_inputs[1]=last_orbit
            date1=date[orbit[first_orbit]:orbit[last_orbit+1]]
            allspec1=allspec[orbit[first_orbit]:orbit[last_orbit+1],:,:] 
            # date1=date[first_orbit:last_orbit-1]
            # allspecextract1=allspecextract[first_orbit:last_orbit-1,:,:] 
            # STOP, 'change eclipse2017 back to orbit'
            allspec1=np.sum(allspec1,axis=1)
            light = np.sum(allspec1,axis=2)
      
            if ploton==True:
                plt.plot(date1, light/max(light),'o')
                plt.xlabel('MJD')
                plt.ylabel('Total Flux')
                plt.show(block=False)
                
            ans = raw_input("Is this correct? (Y/N): ")
            if ans.lower() in ['y','yes']: check=False
            if ploton==True:  plt.close()

    props[1]=eclipse_time(date1, props)
    user_inputs[4]=props[1]
    
    #  We are only interested in scatter within orbits, so correct for flux
    #  between orbits by setting the median of each orbit to the median of
    #  the first orbit
  
    light_corrected=correction(props, date1, light, transit)
  

    # Do a 4-pass sigma cut. 3-5 sigma is ideal. Change n to see how data
    # is affected. A sigma of 3, 4, or 5 could be used, it depends on the
    # data 
    light2=light.copy()
    allspece2=allspec1.copy()
    date2=date1.copy()
    light_corrected2=light_corrected.copy()
    check2=True
    ans2=''

    if check==False:
        light, allspec1, light_corrected, date1 = remove_bad_data(light, allspec1
                                                                  , light_corrected, date1, user_inputs)
    if check==True:
        while check==True:
            light=light2.copy()
            allspec1=allspec2.copy()
            date1=date2.copy()
            light_corrected=light_corrected2.copy()
            
            # This performs the sigma cut and returns input for the fitter: a
            # double array which contains a spectra for each data point
            light, allspec1, light_corrected, date1 = remove_bad_data(light, allspec1
                                                                      , light_corrected, date1
                                                                      , user_inputs, check=check)
            if ploton==True:
                plt.plot(date2, light2,'ro')
                plt.xlabel('MJD')
                plt.ylabel('Total Flux')
                plt.plot(date1, light, 'o',ls='dotted')
                plt.show(block=False)
            ans2=raw_input('This is the new data, with the red points removed. Is this okay? (Y/N): ')
            if ploton==True: plt.close()
            if ans2.lower() in ['y','yes']: check2=False

    results=wl.whitelight2018(props, date1, allspec1, plotting=True, fixtime=True)

    if savedata:

        cols=['Date', 'Spectra']
        processed_data=pd.DataFrame(np.vstack(date1, allspec1),columns=cols)
        processed_data['Visit']=visit+direction
        processed_data['Transit']=transit
        processed_data=processed_data.set_index('Visit')

        sys_p=pd.DataFrame(props, columns=['Properties'])
        sys_p['Visit']=visit+direction
        sys_p=sys_p.set_index('Visit')
 
        try:
            cur=pd.read_csv('./processed_data.csv')
            cur=pd.concat((cur,processed_data))
            cur.to_csv('./processed_data.csv', index=False)
        except IOError:
            processed_data.to_csv('./processed_data.csv', index=False)    
        try:
            curr=pd.read_csv('./system_params.csv')
            curr=pd.concat((curr,sys_p))
            curr.to_csv('./system_params.csv', index=False)
        except IOError:
            sys_p.to_csv('./system_params.csv', index=False)

    return [results, user_inputs]

if __name__==__main__:
    
    if len(sys.argv) < 4:
        sys.exit('Format: run.py [planet] [visit] [direction]')
    visit=sys.argv[1]+'/'+sys.argv[2]
    direction=sys.argv[3]
    transit=False
    if len(sys.argv)==5:
        transit=bool(int(sys.argv[4]))

    best_results, inputs = pre.preprocess_whitelight(visit, direction
                                                     ,transit=transit, savedata=True)

    print "Marg Depth: %d +/- %d" % (best_results[0], best_results[1])
    print "Marg Central Event Time: %d +/- %d" % (best_results[2], best_results[3])
    print "Marg Inclination: %d +/- %d" % (best_results[4], best_results[5])

    inp=pd.DataFrame(inputs, columns='User Inputs')
    inp['Visit']=visit+direction
    inp['Transit']=transit
    inp=inp.set_index('Visit')
    try:
        cur=pd.read_csv('./preprocess_info.csv')
        cur=pd.concat((curr,inp))
        cur.to_csv('./preprocess_info.csv', index=False)
    except IOError:
        inp.to_csv('./preprocess_info.csv', index=False)
 
