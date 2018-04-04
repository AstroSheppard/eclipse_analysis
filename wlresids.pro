PRO wlresids, month, direction, file

IF file eq 'kelt1f1' THEN BEGIN
   xb=0
   start=9-xb
   finish=18-xb
   cc=58-xb
   xc=58
   planet='kelt1'
ENDIF
IF file eq 'kelt1r1' THEN BEGIN
   xb=0
   start=9-xb
   finish=17-xb
   cc=53-xb
   xc=53-xb
   planet='kelt1'
ENDIF
IF file eq 'wasp19a' THEN BEGIN
   start=26-xb ;First exposure of first orbit in whitelight
   finish=55-xb ;last exposure "
   cc=-1-xb ; last exposure used in whitelight: -1 if no extra orbits are ignored
   xb=0 ; First exposure to be included in spectral analysis: 0 if not phase curve
   xc=-1 ; Last "
   planet='wasp19'
ENDIF
IF file eq 'wasp79' THEN BEGIN
   start=12
   finish=24
   cc=-1
   xb=0
   xc=-1
ENDIF
IF file eq 'hatp41' THEN BEGIN
   xb=0
   start=17-xb
   finish=35-xb
   cc=-1
   xc=-1
 ENDIF
 IF file eq 'juner' THEN BEGIN
    xb=0
    start=8-xb
    finish=16-xb
    cc=-1-xb
    xc=-1
 ENDIF
 IF file eq 'junef' THEN BEGIN
    xb=0
    start=9-xb
    finish=17-xb
    cc=-1
    xc=-1
 ENDIF
 IF  file eq 'aug1r' OR file eq 'augusttest2' THEN BEGIN
    xb=0
    start=15-xb
    finish=22-xb
    cc=63-xb
    xc=72
 ENDIF
 IF  file eq 'aug1f' THEN BEGIN
    xb=0
    start=14-xb
    finish=21-xb
    cc=57-xb
    xc=72
 ENDIF
 IF  file eq 'aug2r' THEN BEGIN
    xb=107
    start=113-xb
    finish=118-xb
    cc=-1
    xc=137
 ENDIF
 IF file eq 'aug2f' THEN BEGIN
    xb=100
    start=109-xb
    finish=116-xb
    cc=-1
    xc=139
 ENDIF
IF file eq 'wasp79f' THEN BEGIN
   xb=0
   start=12
   finish=24
   cc=-1
   xc=-1
ENDIF
IF file eq 'aprilf' THEN BEGIN
   xb=0
   start=11
   finish=22
   cc=-1
   xc=-1
ENDIF
IF file eq 'aprilr' THEN BEGIN
   xb=0
   start=11
   finish=21
   cc=-1
   xc=-1
ENDIF
IF file eq 'mayf' THEN BEGIN
   xb=0
   start=11
   finish=22
   cc=-1
   xc=-1
ENDIF
IF file eq 'mayr' THEN BEGIN
   xb=0
   start=11
   finish=21
   cc=-1
   xc=-1
ENDIF



; DEFINE CONSTANTS
  constant = [2.5,20.2,6.67259D-11,2400000.5,86400,7.15D7,6.96D8,1.9D27,1.99D30,5781.6,0.06691666]

; READ IN ALL ZAPPED EXPOSURES

  folder = './zapped2017/' + month + '/' + direction + '/final/*.zap.fits'
  data = FILE_SEARCH(folder)
  nexposure = n_elements(data)
  date=dblarr(nexposure)
  test_expo=MRDFITS(data[0],0,header,/SILENT)
  
; READ IN X AND Y APERTURES
  ap=dblarr(2)
  ; FOR NOW, USE MAX APERTURE
 ; aps, file, ap
 ; x=ap[0]
 ; y=ap[1]
  x=0
  y=0
  xlen = n_elements(test_expo(*,0))+2*x
  ylen = n_elements(test_expo(0,*))+2*y
  expos=dblarr(nexposure,xlen, ylen)

  xmin=-1*x
  xmax=xlen-1-x
  ymin=-1*y
  ymax=ylen-1-y

; RETRIEVE EXPOSURES AND DATES

  FOR img=0, nexposure-1 DO BEGIN
   ; Recover dates
     exposure=MRDFITS(data[img],0,header,/SILENT)
     date[img]=(fxpar(header,'EXPSTART')+fxpar(header,'EXPEND'))/2.
   ;  Cut down image to just what's within aperture
     expo=exposure[xmin:xmax, ymin:ymax]
     expos[img,*,*]=expo
  ENDFOR

  date_order=sort(date)
  date=date[date_order]
  expos=expos[date_order,*,*]
  expos=total(expos,3)
  
  date=date[xb:xc]          
  expos=expos[xb:xc, *]
 
  testbin=expos[*,50:60]

  nexposure=n_elements(date)
;CALCULATE TOTAL FLUX AT EACH TIME
;  start=8 ; starting index change for each month
;  finish=16
  flux = DBLARR(nexposure) 
  flux=TOTAL(expos,2)
  binflux=total(testbin,2)
;  FOR i = 0, n_elements(date)-1 DO BEGIN
;     flux(i) = TOTAL(expos(i,*))
;  ENDFOR
  err = SQRT(flux) 
  binerr=sqrt(binflux)
  binerr = binerr/MEDIAN(binflux[start:finish]) 
 ; err=err * 5.0              
  binflux = binflux/MEDIAN(binflux[start:finish]) ; For june
  errnorm = err/MEDIAN(flux[start:finish]) 
 ; err=err * 5.0              
  fluxnorm = flux/MEDIAN(flux[start:finish]) ; For june


;CALCULATE THE SHIFT IN DELTA_lambda
  sh = DBLARR(nexposure)
 
  FOR i = 0, nexposure-1 DO BEGIN
     sh[i] = CCPEAK(expos[cc,*], expos[i,*],2) 
  ENDFOR

; READ IN SYSTEMATIC MODEL PARAMETERS FROM WHITELIGHT FIT

  parameters=dblarr(24)
  rl=0
  params, file, parameters 
  nModels=n_elements(parameters[*,0])
  radius, file, rl

; CALCULATE HST PHASE AT EACH TIME

  HSTphase = ((date)-(date[start]))/(constant(10))   
  phase2 = FLOOR(HSTphase)
  HSTphase = HSTphase - phase2
  k = WHERE(HSTphase GT 0.6)
 ; kk=WHERE(HSTphase LT -.5)
  IF (k(0) NE -1) THEN HSTphase(k) = HSTphase(k) - 1.0D0
  
; DEFINE ARRAY FOR WHICH TO SAVE WHITELIGHT RESIDUALS

  sys_residuals=dblarr(nModels, nexposure)

;LOOP OVER MODELS

  FOR s=0, nModels -1 DO BEGIN
     
     par=reform(parameters[s,*])
    
;CALCULATE PHASE AT EACH POINT FOR GIVEN MODEL

     phase = ((date)-(par[2]))/par[22]
     phase2 = FLOOR(phase)
     phase = phase - phase2
     a = WHERE(phase GT 0.5)
     IF (a(0) NE -1) THEN phase(a) = phase(a) - 1.0D0
     x2 = phase
  
; CALCULATE ECLIPSE MODEL
     
     b0 = (constant(2)*par(22)*par(22)*86400D0*86400D0/(4*!pi*!pi))^(1D0/3D0) * (par(17)^(1D0/3D0)) * [(sin(x2*2*!pi))^2D0+(cos(par(16))*cos(x2*2*!pi))^2D0]^0.5D0
     plotquery = 0
 
     occultnl, rl, par(18), par(19), par(20), par(21), b0, mulimb0, mulimbf, plotquery

     transit=mulimb0-1d0
     eclipse=1d0 + par(0)*transit/(MAX(transit)-MIN(transit))
 
; CALCULATE SYSTEMATIC MODEL AT EVERY POINT

     systematic_model = (phase*par(3) + 1.0) * (HSTphase*par(4) + HSTphase^2.*par(5) + HSTphase^3.*par(6) + HSTphase^4.*par(7) + HSTphase^5.*par(8) + HSTphase^6.*par(9) + 1.0) * (sh*par(10) + sh^2.*par(11) + sh^3.*par(12) + sh^4.*par(13) + sh^5.*par(14) + sh^6.*par(15) + 1.0)

     w_model = eclipse  * par[1] * systematic_model
     w_cor=fluxnorm/(par[1]*systematic_model)
; CALCULATE RESIDUALS AT EACH POINT

     resids=(fluxnorm-w_model)/par[1]

; SAVE RESIDUALS
     IF s EQ 44 THEN BEGIN
        p=plot(phase, 1- resids, linestyle='', symbol='o', color='purple', position=[.15,.08,.9,.4], yminor=0,ytitle='Model - Data',xtitle='Planetary Phase', ymajor=4, xminor=0)
        p=PLOT(phase, fluxnorm, linestyle='', color='blue',symbol='square', ytitle='Normalized Flux', position=[.15, .4, .9, .95], xshowtext=0, /current, yrange=[.9905,1.002], xminor=0)
        p=PLOT(phase, w_model, symbol='x', linestyle='', /overplot)
      
    ;    t1=text(.5,.60,'White-light Curve', font_size=13, color='blue')
    ;    t2=text(.5,.55,'Best Systematic Model', font_size=13)
    ;   t1=text(.5,.30,'Residuals', font_size=13, color='purple')

;p=plot(wavelength, resid, xtitle='Wavelength [$\mu$m]', ytitle='Residuals', color='green', position=[.15,.08,.9,.3], yminor=0, yrange=[-.003,0.003], xrange=[1.05,1.7],ymajor=5)

;spec=plot(wavelength, spec, ytitle='Electrons', color=555555, /CURRENT,position=[.15,.3,.9,.95], yrange=[2.5D6, 1D7], xshowtext=0, xrange=[1.05,1.7]) 

;a=p.axes

;a[2].hide=1



      
  ;      p=PLOT(phase, binflux, linestyle='', color='blue',symbol='square',xtitle='Planetary Phase', ytitle='Normalized Flux')
  ;      p=PLOT(phase, binflux-.85*resids, linestyle='', color='red',symbol='o',/overplot)
   ;     t1=text(.4,.22,'Spectral Light Curve', font_size=13, color='blue')
   ;     t2=text(.4,.28,'Spectral Light Curve - Residuals', font_size=13, color='red')
  ;      t3=text(.6,.78,'Residuals', font_size=13,color='orchid')
       ; p.close
    
  ;   p=plot(phase, w_cor, symbol='o', linestyle='' , xrange=[-.1,.35], yrange=[.9985, 1.0005])
  ;   p=plot(phase, eclipse, linestyle='-', color='red',/OVERPLOT)
  ;   STOP
   ;  p.close
     ENDIF
     sys_residuals[s,*]=resids
;STOP
 
  ENDFOR
  
; SAVE RESULTS
  first=start
  savename= './WLresiduals/'+file+'.sav'
  SAVE, filename=savename, sys_residuals, expos, sh, hstphase, date, first, cc
  SAVE, filename='./paper/wlcurve/'+file+'.sav', flux, fluxnorm, err, errnorm, date, phase, xc, cc, start, finish, xb
  SAVE, filename='./paper/wlerror/'+planet+'/'+file+'.sav', flux, fluxnorm, err, errnorm, date, phase, xc, cc, start, finish, xb
END

PRO aps, file, ap
  apfile='./apertures2017/'+file+'.sav'
  restore, apfile               ;xbest, ybest, 
  ap[0]=xbest
  ap[1]=ybest
END

PRO params, file, par
  infile='./WLresults/'+file+'.sav'
  RESTORE, infile
  par=sys_params

END

PRO radius, file, rad
  infile='./WLresults/'+file+'.sav'
  RESTORE, infile
  rad=rl

END


; Procedures needed to get eclipse model

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

;+
; NAME:
;           CCPEAK
;
; PURPOSE:
;       Locates the precise location of the peak in the
;       cross-correlation function between two vectors.
;       (Locates LAG_max)
;
; CALLING SEQUENCE:
;
;       LAG_max = CCPEAK(VEC1, VEC2 [, RADIUS ])
;
; INPUTS:
;
;       VEC1, VEC2 - Functions to be cross-correlated

; OPTIONAL INPUTS:
;
;       RADIUS -  How many array elements around the 
;                 nominal peak where polynomial fit 
;                 should be performed.
;
; OUTPUTS:
;
;       LAG_max - Lag at which cross-correlation 
;                 function is maximized
;
; RESTRICTIONS:
;
;       Uses my POLYFIT procedure, not IDL's
;       POLY_FIT
;
;       Uses C_CORRELATE to perform cross-correlation
;
; MODIFICATION HISTORY:
; Written sometime in December 2002 by JohnJohn
; 07 June 2003 JJ - Fixed bug where CF was being indexed out of
; range. Also limited the minimum and maximum lag returned to be
; min(lag) and max(lag), respectively.
; 26 June 2003 JJ - Default radius is now 50 rather than 10
;-

function ccpeak,arr1, arr2, radius, ccf=cf, lag=lag
on_error,2		;Return to caller if an error occurs
n = n_elements(arr1)
if n_elements(radius) eq 0 then radius = 50
lag = fillarr(1,-radius,radius)
cf = c_correlate(arr1,arr2,lag)
dum = max(cf, ind)

srad = 3
sublag = lag[(ind-srad) > 0:(ind+srad) < (2*radius)]
subcf  = cf[(ind-srad) > 0:(ind+srad) < (2*radius)]
a = polyfit(sublag,subcf,2)
maxlag = - a[1]/(2.*a[2])
nlag = n_elements(lag)
if maxlag lt lag[0] then maxlag = lag[0]
if maxlag gt lag[nlag-1] then maxlag = lag[nlag-1]
return, maxlag
end

;+
; NAME: 
;           FAN
;
; PURPOSE:
;           Take the outer product of the input ARRAY and a
;           UNIT_VECTOR to "fan out" a 1D vector into an array 
;           comprised of the vector repeated row-wise NFAN times.
;           Useful for array-wise mathematics (Look Ma, no FOR loops!)
;
; CALLING SEQUENCE:
;           result = fan(array [,nfan, /transpose])
;
; INPUTS:
;           ARRAY - 1D array, input vector
;           NFAN  - number of times to repeat the input vector,
;                   default is N_ELEMENTS(ARRAY)
;
; KEYWORD PARAMETERS:
;                     
;           TRANSPOSE - Repeat the input vector column-wise
;
; OUTPUTS:
;           A 2D array with N_ELEMENTS(ARRAY) columns and NFAN
;           rows.
;
; EXAMPLE:
;           Fan a FINDGEN of 3 elements, twice.
;
;           IDL> a = findgen(3)
;           IDL> print,fan(a,2)
;                 0.00000      1.00000      2.00000
;                 0.00000      1.00000      2.00000
;
; MODIFICATION HISTORY:
;           Created sometime in ought-2 by JohnJohn
; 06 Dec 2002 JohnJohn- Added some error handling at the beginning
;-
function fan,array,nfan,transpose=transpose
  on_error,2                    ;if broke then return to sender
  if n_params() lt 1 then begin 
     message,'Syntax: f = fan(array [,nfan, /transpose])',/info
     return,-1
  endif

  if n_elements(nfan) eq 0 then nfan = n_elements(array)
  unit_vector = replicate(1d,nfan) ;dblarr(nfan)+1.
  if keyword_set(transpose) then new = array##unit_vector $
  else new = unit_vector##array
  return,new
end

;+
; NAME: 
;       FILLARR 
;
;
; PURPOSE:
;       This function generates an array from MIN to MAX with
;       step size DEL. If an integer number of steps cannot be
;       fit between MIN and MAX, then MAX will be adjusted to
;       be as close as the specified maximum as possible.
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;       f = fillarr(n, min, max [,fan=, transfan=, /double])
;
;
; INPUTS:
;       DEL:  The desired step size
;       MIN:  The value of the first array element in F
;       MAX:  The value of the last array element in F if
;             (MAX-MIN)/DEL is an integer. Adjusted otherwise.
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
;       FANNNED:    Number of times the array is to be repeated.
;                   The final dimensions of F  will be 
;                   fix((MAX-MIN)/DEL) + 1 columns by FANNED ows.
;
;       /TRANSPOSE  Final dimensions of F wil be FAN columns by 
;                   fix((MAX-MIN)/DEL) + 1 rows if FAN is specified. 
;
; OUTPUTS:
;
;       F:    Final array. If input parameters are double precision,
;             then F will be double as well. F is float otherwise.
;
; RESTRICTIONS:
;
;       You'll need FAN.PRO to use the fan= keyword. 
;       http://astron.berkeley.edu/~johnjohn/idl.html#FAN
;
; EXAMPLE:
;
;         For an array that runs from 2 to 5 in steps of .7
;
;         IDL> f = fillarr(.7,2,5)
;         IDL> print, f
;            2.00000      2.70000      3.40000     4.10000    4.80000
;         
; MODIFICATION HISTORY:
; Written by John "JohnJohn" Johnson 21-Feb-2002
; 22-Feb-2002 JohnJohn- Fixed precision bug
; 23-Feb-2002 JohnJohn- Calculations performed in double precision. 
;                       Output in double precision if input is 
;                       double.
; 01-Mar-2002 JohnJohn- Much props to Tim Robishaw (Tdogg) for helping
;                       me understand machine precision and finally fixing
;                       the precision bug.
; 23 Apr 2002 JohnJohn- Modified the /FAN operation to match my new
;                       FAN procedure. Default direction of the
;                       fanning process is in the column direction,
;                       i.e. a 5-element array with FAN=2 will yeild a
;                       5x2 array rather than the other way around.
; 06 Dec 2002 JohnJohn- Modified the /FAN operation again to run using
;                       the actuall FAN proceedure which is faster
;                       than doing two separate operations for fanning
;                       and taking the transpose. duh.
; 14 Apr 2005 JohnJohn- Fixed bug where if n_params() eq 2, then MIN
;                       was being modified. Input variable now
;                       protected by renaming MIN as MININ.
;-
function fillarr,del,minin,max,fanned=fanned,transpose=transpose
;DEAL WITH HUMANS
on_error,2		;Return to caller if an error occurs
if n_params() lt 2 then message,'INCORRECT NUMBER OF INPUTS. Syntax: f = fillarr(del,min,max)',/ioerror

if n_params() eq 2 then begin
    max = minin[1]
    min = minin[0]
endif else min = minin

if max lt min then message,'MIN must be less than MAX',/ioerror
if del eq 0 then message,'DEL cannot equal 0',/ioerror

;if all of the input parameters are double, the return the answer in
;double precision.
doub = (size(del,/type) eq 5) and (size(min,/type) eq 5) and (size(max,/type) eq 5) or keyword_set(double)
del = double(del)
min = double(min)
max = double(max)
;ARG will go into A later. These are the only real calculations performed.
arg = (max-min)/del
;test for and correct rounding errors
rnd = round(arg)
eps = (machar(/double)).eps
if abs(rnd-arg) lt rnd*eps then arg = rnd else arg = fix(arg,type=3)

a = dindgen(arg+1)*del+min      ;can you believe there's all this code just to do this?

if n_elements(fanned) ne 0 then begin
    nfan = fanned
    if keyword_set(transpose) then a = fan(a,nfan, /transpose) $
    else a = fan(a,nfan)
endif

if not doub then a = float(a)

return,a 
end

function polyfit,t,y,deg,yfit,yfit1,covariance=cov,weight=w
;on_error,2
n = n_elements(t)
pow = indgen(deg+1)
powarr = fan(pow,n,/trans)
x =  fan(double(t),deg+1)
xarr = x^powarr
xarrt = transpose(xarr)
if keyword_set(w) then xarr = fan(double(w),deg+1)*x^powarr 

alpha = xarr##xarrt
beta = xarr##(double(y))
cov = invert(alpha)
a = cov##beta
if n_params() eq 4 then yfit = poly(t,a)
if n_params() eq 5 then yfit1 = poly(t,a) ;;Legacy code, kept for compatibility
return,a
end
