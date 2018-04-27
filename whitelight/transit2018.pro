PRO transit2018, month, direction, X=x, Y=y, results, user_inputs,PLOTON=plotting, CHECK=checks, INP_FILE=inp,SAVEFILE=name, SAVEDATA=savedata

; INPUTS
 ; Month/Direction allows the user to access the data
 ; x, y, allow the user to input a certain aperture
  ; x and y are both 1-D 3 elements arrays: [aperture, min pixel, max
  ; pixel]
  ; residuals: Returns the standard deviation of the residuals for the
  ; best model
 ; plotting: set to "on" to see plots needed to check data
 ; checks: set to "on" to manually reduce data
 ; If checks is set to on, "user_inputs" will return the inputs
  ; that the user used: [first orbit, last orbit, sigma cut factor, 
  ; number of passes, center eclipse time]. If checks is set to off, then
  ; the user_inputs array will be used as inputs (easier to automate)
 
; This does it all from scratch. If data is
;already zapped, leave this commented out to save time
; Add keyword for zapping
;wasp18_zapping, month 
set_plot, 'x'
folder = './zapped2017/' + month + '/' + direction + '/final/*.zap.fits'
data = FILE_SEARCH(folder)
;; Locate all fo the fits images to extract - these are all of the
;; zapped image files: they only contain source and have removed
;; cosmic rays

; can ignore this now. If no x provided, just use xlen and ylen
;IF n_elements(x) EQ 0 THEN BEGIN
;   IF n_elements(inp) EQ 0 THEN BEGIN

;      x=dblarr(3)               ;[apt,min,max]
;      y=dblarr(3)
;      APERTURE, data, nexposure, x, y
;   ENDIF ELSE BEGIN
;      save_folder='/dosquadequis1/ksheppard/Desktop/2016_research/sav_files/'
;      RESTORE, save_folder+inp
;      x=xbest
;      y=ybest
;   ENDELSE
;ENDIF
 ;check x, y, inputs

nexposure = n_elements(data)
print, 'There are ', nexposure, '     exposures in this visit'

date=dblarr(nexposure)
time=dblarr(nexposure) ; comment out

test_expo=MRDFITS(data[0],0,header,/SILENT)
xlen = n_elements(test_expo(*,0))+2*x
ylen = n_elements(test_expo(0,*))+2*y
allspecextract=dblarr(nexposure,xlen, ylen)

xmin=-1*x
xmax=xlen-1-x
ymin=-1*y
ymax=ylen-1-y

FOR img=0, nexposure-1 DO BEGIN

   ; Recover dates
   exposure=MRDFITS(data[img],0,header,/SILENT)
   date[img]=(fxpar(header,'EXPSTART')+fxpar(header,'EXPEND'))/2.
   time[img]=fxpar(header,'EXPTIME')
; Cut down image to just what's within aperture
   expo=exposure[xmin:xmax, ymin:ymax]
 ;  zapping, x[0], y[0], expo, zapped, 'off'
   allspecextract[img,*,*]=expo
ENDFOR
orbit = dblarr(1)
; Change this so input is read in from program
inputs, './'+month+'/inputs.dat', props

;input=[0.05917,0,1.500983,3901.91314,0.9414523,1] 
folder='./' 



 ; ops, FILE=month+".ps"


date_order=sort(date)
date=date[date_order]
allspecextract=allspecextract[date_order,*,*]

  ; Classify the data by each HST orbit. Returns array (orbit) 
  ; which contains the indeces for the start of each orbit

ORBITS, date, orbit
print
print, "Number of total orbits: ", n_elements(orbit)

  ; Choose which orbits to include in the eclipse fitting. 1-2 on either
  ; side of the eclipse is recommended

check = 1
ans=''

IF ~KEYWORD_SET(checks) THEN BEGIN
   IF n_elements(inp) EQ 1 AND n_elements(user_inputs) EQ 0 THEN BEGIN
      RESTORE, folder+'apertures2017/'+inp
      user_inputs=u_inputs
   ENDIF
   IF n_elements(user_inputs) EQ 0 AND n_elements(inp) EQ 0 THEN BEGIN
      STOP, 'You must either manually pick the inputs, load them from another program, or upload them from a saved file. Provide a file or set /CHECK and /PLOT.'
   ENDIF

   first_orbit=user_inputs[0]
   last_orbit=user_inputs[1]
   date1=date[orbit[first_orbit]:orbit[last_orbit]-1]
   allspecextract1=allspecextract[orbit[first_orbit]:orbit[last_orbit]-1,*,*] 
   allspecextract1=total(allspecextract1,3) ;spectra for each exposure
   light = total(allspecextract1, 2) ; total light for each exposure
  

ENDIF


IF KEYWORD_SET(checks) THEN BEGIN
   user_inputs=dblarr(5)
   WHILE check DO BEGIN
      IF KEYWORD_SET(plotting) THEN plot=PLOT(date, total(total(allspecextract,3),2),  xtitle='MJD', ytitle='Total Flux',symbol='o',linestyle='')
      READ, first, PROMPT="Enter the first orbit to include (starting from 0): "
      first_orbit=FIX(first)
      user_inputs[0]=first_orbit
      READ, last, PROMPT="Enter the last orbit to include (starting form 0): "
      last_orbit=FIX(last)+1
      IF KEYWORD_SET(plotting) THEN plot.Close
      user_inputs[1]=last_orbit
      date1=date[orbit[first_orbit]:orbit[last_orbit]-1]
      allspecextract1=allspecextract[orbit[first_orbit]:orbit[last_orbit]-1,*,*] 
     ; date1=date[first_orbit:last_orbit-1]
     ; allspecextract1=allspecextract[first_orbit:last_orbit-1,*,*] 
     ; STOP, 'change eclipse2017 back to orbit'
      allspecextract1=total(allspecextract1,3)
      light = total(allspecextract1, 2)
      
      IF KEYWORD_SET(plotting) THEN BEGIN
         plot=PLOT(date1, light/max(light),  xtitle='MJD', ytitle='Total Flux',symbol='o',linestyle='')
      ENDIF
      READ, ans, PROMPT="Is this correct? (Y/N): "
      IF (ans EQ 'Y') THEN check=0
      IF KEYWORD_SET(plotting) THEN  plot.Close
   ENDWHILE
ENDIF

eclipse_time, date1, props
print, date1
stop, props[1]
; We are only interested in scatter within orbits, so correct for flux
; between orbits by setting the median of each orbit to the median of
; the first orbit
                            ; check x, y, user_inputs, input
 ; Note: Only works for 3-5 orbits currently
  
;CORRECTION, first_orbit, last_orbit, orbit, light, light_corrected
CORRECTION, props, date1, light, light_corrected


; Do a 4-pass sigma cut. 3-5 sigma is ideal. Change n to see how data
; is affected. A sigma of 3, 4, or 5 could be used, it depends on the
; data 
light2=light
allspecextract2=allspecextract1
date2=date1
light_corrected2=light_corrected
check2=1
ans2=''

IF ~KEYWORD_SET(checks) THEN REMOVE_BAD_DATA, light, allspecextract1, light_corrected, date1, user_inputs
IF KEYWORD_SET(checks) THEN BEGIN
   WHILE check2 DO BEGIN
      light=light2
      allspecextract1=allspecextract2
      date1=date2
      light_corrected=light_corrected2

     
      
; This performs the sigma cut and returns input for the fitter: a
; double array which contains a spectra for each data point

      REMOVE_BAD_DATA, light, allspecextract1, light_corrected, date1, user_inputs, /CHECK
      IF KEYWORD_SET(plotting) THEN BEGIN
         plot2=plot(date2, light2,  xtitle='MJD', ytitle='Total Flux',symbol='o', color='r', linestyle='')
         plot2=plot(date1, light,  xtitle='MJD', ytitle='Total Flux',symbol='o',linestyle=1, /OVERPLOT)
      ENDIF
      READ, ans2, PROMPT='This is the new data, with the red points removed. Is this okay? (Y/N): '
      IF KEYWORD_SET(plotting) THEN plot2.Close
      IF (ans2 EQ 'Y') THEN check2=0
   ENDWHILE
ENDIF
                            ; check x, y, user_inputs, input
  ; Determine the size of the first orbit
orbit2=dblarr(1)
ORBITS, date1, orbit2
first_orbit_size=orbit2[1]

error=stddev(light_corrected)/median(light_corrected[0:first_orbit_size-1])


; HERE
eclipse_time, date1, props
user_inputs[4]=props[1]
results=dblarr(9)
;IF KEYWORD_SET(name) THEN save_file = folder + 'WL/'+ name
IF ~KEYWORD_SET(name) THEN name= STRMID(month,0,3) + STRMID(direction,0,1)

;STOP                           ; check x, y, user_inputs, input

  
;IF KEYWORD_SET(name) THEN BEGIN
;   whitelight_eclipse, props, date1, allspecextract1,first_orbit_size,results, SAVEFILE=save_file
;ENDIF ELSE BEGIN
;save, filename=folder+'august.sav', allspecextract1, date1
;whitelight_eclipse_err, props, date1, allspecextract1, first_orbit_size,error, /FIXTIME,  saveerr=month+direction

; Have to toggle FIXTIME based on whether it can eb constrained by data
whitelight_eclipse2017, props, date1, allspecextract1,first_orbit_size,results, /PLOTTING, /FIXTIME;, SAVEFILE=name
;ENDELSE



IF KEYWORD_SET(savedata) THEN BEGIN
    save_folder=folder+'eclipse_inputs/'
    savefile=save_folder+name+'.sav'
    SAVE, filename=savefile, props, date1, allspecextract1, first_orbit_size, results
ENDIF
;cps

END






; Procedure to remove columns from a 2D array

FUNCTION RemoveCols, array, cols

   ; array -- A 2D array from which rows will be removed.
   ; rows -- A vector of row indices to remove from array.

   ; Need both positional parameters.
   IF N_Params() NE 2 THEN BEGIN
       Print, "Usage: 'RemoveRows, array, rowsToRemove'"
       RETURN, -1
    ENDIF
    
    ; The array must be 2D.
    ndims = Size(array, /N_DIMENSIONS)
    IF ndims NE 2 THEN BEGIN
        void = Dialog_Message('Array must be 2D.')
        Print, "Usage: 'RemoveRows, array, rowsToRemove'"
        RETURN, -1
    ENDIF
    
    ; The rows must be a vector.
    IF Size(cols, /N_DIMENSIONS) EQ 0 THEN cols = [cols]
    
    ; Find the dimensions of the array.
    dims = Size(array, /DIMENSIONS)

    ; Return the shortened array.
    RETURN, array[Where(~Histogram(cols, MIN=0, MAX=dims[0]-1), /NULL),*]
        
END


; Procedure to organize light curve data by HST orbit

PRO ORBITS, date, orbit
  FOR i=0, n_elements(date)-2 DO BEGIN
     t=date[i+1]-date[i]
     IF t*86400 GT 1200 THEN orbit=[orbit, i+1] ; 1800s is ~ half an HST orbit
  ENDFOR
  orbit=[orbit,n_elements(date)]
END



;Procedure to remove "bad" data from light curve

PRO REMOVE_BAD_DATA, light_curve, model_input, light_corrected, date1, CHECK=checks, user_inputs
; Plot the corrected light curve with 3-4-5 sigma marked
  med= median(light_corrected)
  sigma = SQRT(TOTAL((light_corrected-med)^2)/(2*n_elements(light_corrected)))
 
  medi=dblarr(n_elements(date1))
  sig3=dblarr(n_elements(date1))
  sig4=dblarr(n_elements(date1))
  sig5=dblarr(n_elements(date1))
  medi[*]=med
  sig3[*]=med+3*sigma
  sig4[*]=med+4*sigma
  sig5[*]=med+5*sigma

  IF ~KEYWORD_SET(checks) THEN BEGIN
     nPasses=user_inputs[3]
     sigma_cut_factor=user_inputs[2]
  ENDIF

  IF KEYWORD_SET(checks) THEN BEGIN
     plot=plot(date1, light_corrected,  xtitle='MJD', ytitle='Total Flux',symbol='o',linestyle=1)
     plot=plot(date1, sig5,color=113454322,/OVERPLOT, NAME= '5-sigma')
     plot1=plot(date1, sig4,color=1224332,/OVERPLOT, NAME= '4-sigma')
     plot2=plot(date1, sig3,color=13312532,/OVERPLOT, NAME= '3 sigma')
     plot3=plot(date1, medi,/OVERPLOT, NAME='Median')
     leg=LEGEND(TARGET=[plot,plot1,plot2,plot3],/AUTO_TEXT_COLOR)
 
; Ask the user for a sigma cut factor and the number of passes

     READ, cut, PROMPT="Enter the sigma-cut factor (3-5 recommended): "
     sigma_cut_factor=DOUBLE(cut)
     user_inputs[2]=sigma_cut_factor
     READ, pass, PROMPT="Enter the number of passes for the sigma-cut: "
     nPasses=FIX(pass)
     user_inputs[3]=nPasses
  
     plot.Close
  ENDIF
  
; Cut out the "bad" data

  FOR j=0, nPasses DO BEGIN
     med= median(light_corrected)
     sigma = SQRT(TOTAL((light_corrected-med)^2)/(2*n_elements(light_corrected)))
     dif= ABS(light_corrected-med)
     index=WHERE(dif GT sigma_cut_factor*sigma)
     IF index[0] NE -1 THEN BEGIN
        ; Remove bad data points
        REMOVE, index, light_curve, date1, light_corrected
        model_input=RemoveCols(model_input,index)
     ENDIF
  ENDFOR
END


PRO APERTURE, data, nexposure, x, y
  exposure = MRDFITS(data(0),0,/silent) ; use first frame to get length and height info

  xlen = n_elements(exposure(*,0))
  ylen = n_elements(exposure(0,*))

  exposure_xcen = DBLARR(nexposure,ylen)
  exposure_ycen = DBLARR(nexposure,xlen)
  scan_width = DBLARR(nexposure)
  scan_center = DBLARR(nexposure)
  scan_height = DBLARR(nexposure)
  center=DBLARR(nexposure, 2)
  maxi=dblarr(nexposure)

  FOR img=0, nexposure-1 DO BEGIN
     
     exposure=MRDFITS(data[img],0,header,/silent)
     xlen = n_elements(exposure(*,0)) 
     ylen = n_elements(exposure(0,*))

; x center is just the mipoint in the x-direction of the spectra. Set
; xcen array to 2-D array: 1D for each exposure, 2D for each y value
; at the center x for that exposure
     exposure_xcen(img,*) = exposure((xlen/2.),*) 

     max_count=max(exposure)
     maxi(img)=max_count
;calculate exposure ycen
     scan = WHERE(exposure_xcen(img,*) GT max_count/10.) ; Find y- coordinates of pixels where the (center x, y) pixel value is greater than 50 (bright)

     scan_width(img) = n_elements(scan) ; width of source
     scan_center(img) = scan(0)+(scan_width(img)/2.) ; midpoint between first and last significant y-values is y-center
     exposure_ycen(img,*) = exposure(*,scan_center(img)) ; set another 2D array (1st D images) with the x-value for each x-coordiante with the center y-coordinate

     height = WHERE(exposure_ycen(img,*) GT max_count/10.) ; At center y, scan upward until there is no longer a significant signal. This is the height of the source
     scan_height(img) = n_elements(height)
     height_cen=height[0]+scan_height[img]/2.


     exposure_xcen[img,*]=exposure[height_cen,*]
     scan2=where(exposure_xcen[img,*] GT max_count/10.)
     scan_width[img]= n_elements(scan2)
     scan_center[img]=scan2[0]+scan_width[img]/2.

     center[img,0]=height_cen
     center[img,1]=scan_center[img]

  ENDFOR


  xapt=median(scan_height/2)
  yapt=median(scan_width/2)
  xmin=median(center[*,0]-xapt)
  xmax=median(center[*,0]+xapt)
  ymin=median(center[*,1]-yapt)
  ymax=median(center[*,1]+yapt)

  x[0]=xapt
  x[1]=xmin
  x[2]=xmax
  y[0]=yapt
  y[1]=ymin
  y[2]=ymax
END

PRO inputs, data, props

  OPENR, lun, data, /GET_LUN
  ncols=3
  nrows=FILE_LINES(data)
  label=strarr(nrows)
  units=strarr(nrows)
  data_arr=dblarr(nrows)
  READCOL, data, F='A,A,D',label,units,data_arr
  CLOSE, lun
  FREE_LUN, lun
 
;Rj-m, Rsolar-m,deg-radian,days-secs,AU-m
  conversions=[6.9911e7,6.957e8,!pi/180.,86400,1.49598e11]
  G=6.67408e-11                 ; SI units
  depth=data_arr[0]*conversions[0]/(data_arr[1]*conversions[1])*(data_arr[2]/data_arr[3])^.5
  inc=data_arr[5]*conversions[2]
  epoch=data_arr[6]-2400000.5   ;Last observed eclipse
  period=data_arr[4]
  a_R=data_arr[7]*conversions[4]/(data_arr[1]*conversions[1])
  MpMsR=(a_R^3)*4*!PI*!PI/(G*period*period*conversions[3]*conversions[3])
  
  rl = data_arr[0]*conversions[0]/(data_arr[1]*conversions[1]) ;Rp/Rs
  
  props=dblarr(7)
;  props[0]=depth ;actually sqrt(depth)
  props[0]=rl
  props[1]=epoch
  props[2]=inc
  props[3]=MpMsR
  props[4]=period
  props[5]=1
  props[6]=depth*depth/3 ; actually depth
  
 ; SAVE, FILENAME='properties.sav', props

END

; Program to determine the expected eclipse time
PRO eclipse_time, date, properties
  ; Inputs
  ; date: 1D array of the date of each exposure (MJD)
  ; properties: 1D array containing the last observed eclipse 
  ; and the period. (MJD, days)
  ; NOTE: If transit has been observed but eclipse hasn't, then 
  ; use last observed transit and add period/2
  time=properties[1]
  period=properties[4]
  WHILE (time LT date[0]) DO BEGIN
     time=time+period
  ENDWHILE
  properties[1]=double(time)
END
  
PRO correction, inputs, date, flux, corrected
  t0=inputs[1]
  depth=inputs[6]
  MsMpR=inputs[3]
  inclin=inputs[2]
  Per = inputs[4] 
  rl=inputs[0]
  JD = 2400000.5D0              
  Gr = 6.67259D-11


  phase = (date-t0)/(Per) 
  phase2 = FLOOR(phase)
  phase = phase - phase2
  a = WHERE(phase GT 0.5)
  IF (a(0) NE -1) THEN phase(a) = phase(a)-1.0D0


   b0 = (Gr*Per*Per*86400D0*86400D0/(4*!pi*!pi))^(1D0/3D0) * (MsMpR^(1D0/3D0)) * [(sin(phase*2*!pi))^2D0 + (cos(inclin)*cos(phase*2*!pi))^(2D0)]^(0.5D0)

   plotquery = 0
   
   occultnl,rl,0,0,0,0,b0,mulimb0,mulimbf,plotquery,_extra=e

   transit=mulimb0-1d0
   eclipse=(1d0 + depth*transit/(MAX(transit)-MIN(transit)))
   corrected=flux/eclipse
   
END

; Procedure to nomalize eclipse so that points can be eliminated based
; on scatter. Possibly problematic for orbits in egress or ingress

;PRO CORRECTION, first_orbit, last_orbit, orbit_array, light_curve, light_corrected
;  med=dblarr(last_orbit-first_orbit) ; Median of each orbit
;  n=dblarr(last_orbit-first_orbit) ; # data points per orbit
;  i=0
;  FOR j=first_orbit, last_orbit-1 DO BEGIN
;     orb=light_curve[orbit_array[j]-orbit_array[first_orbit]:orbit_array[j+1]-orbit_array[first_orbit]-1]
;     med[i]=median(orb)
;     n[i]=n_elements(orb)
;     i+=1
;  ENDFOR

;  light_correction = dblarr(n[0])
;  FOR i=1, n_elements(n)-1 DO BEGIN
;     orb=dblarr(n[i])+med[0]-med[i]
;     light_correction=[light_correction, orb]
;  ENDFOR

;  light_corrected=light_curve+light_correction
;END



pro occultnl,rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,plotquery,_extra=e
; Please cite Mandel & Agol (2002) if making use of this routine.
timing=systime(1)
; This routine uses the results for a uniform source to
; compute the lightcurve for a limb-darkened source
; (5-1-02 notes)
;Input:
;  rl        radius of the lens   in units of the source radius
;  c1-c4     limb-darkening coefficients
;  b0        impact parameter normalized to source radius
;  plotquery =1 to plot magnification,  =0 for no plots
;  _extra=e  plot parameters
;Output:
; mulimb0 limb-darkened magnification
; mulimbf lightcurves for each component
; 
; First, make grid in radius:
; Call magnification of uniform source:
occultuniform,b0,rl,mulimb0
bt0=b0
fac=max(abs(mulimb0-1.d0))
;print,rl
omega=4.d0*((1.d0-c1-c2-c3-c4)/4.d0+c1/5.d0+c2/6.d0+c3/7.d0+c4/8.d0)
nb=n_elements(b0)
indx=where(mulimb0 ne 1.d0)
mulimb=mulimb0(indx)
mulimbf=dblarr(nb,5)
mulimbf(*,0)=mulimbf(*,0)+1.d0
mulimbf(*,1)=mulimbf(*,1)+0.8d0
mulimbf(*,2)=mulimbf(*,2)+2.d0/3.d0
mulimbf(*,3)=mulimbf(*,3)+4.d0/7.d0
mulimbf(*,4)=mulimbf(*,4)+0.5d0
nr=2
dmumax=1.d0
;while (dmumax gt fac*1.d-10 and nr le 16) do begin
while (dmumax gt fac*1.d-3) do begin
;while (dmumax gt 1.d-6 and nr le 4) do begin
  mulimbp=mulimb
  nr=nr*2
  dt=0.5d0*!pi/double(nr)
  t=dt*dindgen(nr+1)
  th=t+0.5d0*dt
  r=sin(t)
  sig=sqrt(cos(th(nr-1)))
  mulimbhalf =sig^3*mulimb0(indx)/(1.d0-r(nr-1))
  mulimb1    =sig^4*mulimb0(indx)/(1.d0-r(nr-1))
  mulimb3half=sig^5*mulimb0(indx)/(1.d0-r(nr-1))
  mulimb2    =sig^6*mulimb0(indx)/(1.d0-r(nr-1))
  for i=1,nr-1 do begin
; Calculate uniform magnification at intermediate radii:
    occultuniform,b0(indx)/r(i),rl/r(i),mu
; Equation (29):
    sig1=sqrt(cos(th(i-1)))
    sig2=sqrt(cos(th(i)))
    mulimbhalf =mulimbhalf +r(i)^2*mu*(sig1^3/(r(i)-r(i-1))-sig2^3/(r(i+1)-r(i)))
    mulimb1    =mulimb1    +r(i)^2*mu*(sig1^4/(r(i)-r(i-1))-sig2^4/(r(i+1)-r(i)))
    mulimb3half=mulimb3half+r(i)^2*mu*(sig1^5/(r(i)-r(i-1))-sig2^5/(r(i+1)-r(i)))
    mulimb2    =mulimb2    +r(i)^2*mu*(sig1^6/(r(i)-r(i-1))-sig2^6/(r(i+1)-r(i)))
  endfor
  mulimb=((1.d0-c1-c2-c3-c4)*mulimb0(indx)+c1*mulimbhalf*dt+c2*mulimb1*dt+$
           c3*mulimb3half*dt+c4*mulimb2*dt)/omega
  ix1=where(mulimb+mulimbp ne 0.d0)
  dmumax=max(abs(mulimb(ix1)-mulimbp(ix1))/(mulimb(ix1)+mulimbp(ix1)))
;  print,'Difference ',dmumax,' nr ',nr
endwhile
mulimbf(indx,0)=mulimb0(indx)
mulimbf(indx,1)=mulimbhalf*dt
mulimbf(indx,2)=mulimb1*dt
mulimbf(indx,3)=mulimb3half*dt
mulimbf(indx,4)=mulimb2*dt
mulimb0(indx)=mulimb
if(plotquery eq 1) then plot,bt0,mulimb0,_extra=e
if(plotquery eq 1) then oplot,bt0,mulimbf(*,0),linestyle=2
b0=bt0
;print,'Time ',systime(1)-timing
return
end

pro occultuniform,b0,w,muo1
if(abs(w-0.5d0) lt 1.d-3) then w=0.5d0
; This routine computes the lightcurve for occultation
; of a uniform source without microlensing  (Mandel & Agol 2002).
;Input:
;
; rs   radius of the source (set to unity)
; b0   impact parameter in units of rs
; w    occulting star size in units of rs
;
;Output:
; muo1 fraction of flux at each b0 for a uniform source
;
; Now, compute pure occultation curve:
nb=n_elements(b0)
muo1=dblarr(nb)
for i=0,nb-1 do begin
; substitute z=b0(i) to shorten expressions
z=b0(i)
; the source is unocculted:
; Table 3, I.
if(z ge 1.d0+w) then begin
  muo1(i)=1.d0
  goto,next
endif
; the  source is completely occulted:
; Table 3, II.
if(w ge 1.d0 and z le w-1.d0) then begin
  muo1(i)=0.d0
  goto,next
endif
; the source is partly occulted and the occulting object crosses the limb:
; Equation (26):
if(z ge abs(1.d0-w) and z le 1.d0+w) then begin
  kap1=acos(min([(1.d0-w^2+z^2)/2.d0/z,1.d0]))
  kap0=acos(min([(w^2+z^2-1.d0)/2.d0/w/z,1.d0]))
  lambdae=w^2*kap0+kap1
  lambdae=(lambdae-0.5d0*sqrt(max([4.d0*z^2-(1.d0+z^2-w^2)^2,0.d0])))/!pi
  muo1(i)=1.d0-lambdae
endif
; the occulting object transits the source star (but doesn't
; completely cover it):
if(z le 1.d0-w) then muo1(i)=1.d0-w^2
next:
endfor
;muo1=1.d0-lambdae
return
end


