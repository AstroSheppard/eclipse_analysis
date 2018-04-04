PRO fullzap, month


; don't forget to compile data_zapping.pro
;;;;;; Cosmic ray zap ;;;;;;;;;;
 ; data_zapping, month

;;;;;;;;;; Bad pixel zap ;;;;;;;

; read in zapped data
  folder='./zapped2017/'+month+'/forward/cr/*.fits'
  data=FILE_SEARCH(folder)
  n=n_elements(data)


  x=dblarr(3)                   ;[apt,min,max]
  y=dblarr(3)
  xl=x
  yl=y
  IF n GT 1 THEN BEGIN
;;;;;;; Forward ;;;;;;;;;;;
;Find base aperture 
     aperture, data, n, x, y


; Choose largest tested aperture
     xl=x+4
     xl[1]=x[1]-4
     yl=y+4
     yl[1]=y[1]-4

; Go through each exposure and zap bad pixels

     FOR img=0, n-1 DO BEGIN
        exposure=MRDFITS(data[img],0,/SILENT)
        expo=exposure[xl[1]:xl[2]-1, yl[1]:yl[2]-1]
        header=headfits(data[img])
           
      ;  IF img eq 4 then begin
      ;     STOP
           pixel_zapping, xl[0], yl[0], expo,  zapped, 'off'
                                ; save zapped images in folder named by dataset
           filename = './zapped2017/'+ month + '/forward/final/' + STRING(img,format='(I03)') +'.zap.fits'
           MWRFITS, zapped, filename, header, /CREATE
           print, 'image', img
      ;  ENDIF
     ENDFOR
  ENDIF

; Now, do same for reverse scan
  folder='./zapped2017/'+month+'/reverse/cr/*.fits'
  data=FILE_SEARCH(folder)
  n=n_elements(data)

  x=dblarr(3)                   ;[apt,min,max]
  y=dblarr(3)
  IF n GT 1 THEN BEGIN
     
     aperture, data, n, x, y

     xl=x+4
     xl[1]=x[1]-4
     yl=y+4
     yl[1]=y[1]-4

     FOR img=0, n-1 DO BEGIN
        exposure=MRDFITS(data[img],0,/SILENT)
        expo=exposure[xl[1]:xl[2]-1, yl[1]:yl[2]-1]
        header=headfits(data[img])
        
        pixel_zapping, xl[0], yl[0], expo,  zapped, 'off'
        filename = './zapped2017/'+ month + '/reverse/final/' + STRING(img,format='(I03)') +'.zap.fits'
        MWRFITS, zapped, filename, header, /CREATE
       
     ENDFOR
  ENDIF
  
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
  date=dblarr(nexposure) ; test 
  qual=strarr(nexposure)
  FOR img=0, nexposure-1 DO BEGIN
     
     exposure=MRDFITS(data[img],0,header,/silent)
     xlen = n_elements(exposure(*,0)) 
     print, xlen, img
     ylen = n_elements(exposure(0,*))
; tests
     date[img]=(fxpar(header,'EXPSTART')+fxpar(header,'EXPEND'))/2.
     qual[img]=fxpar(header,'QUALITY')
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

;  xmin=center[*,0]-scan_height/2
;  xmax=center[*,0]+scan_height/2
;  d=sort(date)
;  date=date[d]
;  xmin=xmin[d]
;  xmax=xmax[d]
;  plot=PLOT(date[9:60], xmin[9:56], color='red',symbol='x',linestyle='')
;  plot=PLOT(date[9:60], xmax[9:56], symbol='x',linestyle='',color='red',/OVERPLOT)
;  plot=PLOT(date[9:60], center[9:56,0],symbol='x',linestyle='',color='green',/OVERPLOT)
;  STOP
;  plot.Close
;  STOP
 ; histplot.Close
 ; hist=HISTOGRAM(center, LOCATIONS=xbin, nbins=10)
 ; histplot=PLOT(xbin, hist, xrange=[0,max(xbin)])
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


PRO pixel_zapping, xapt, yapt, allspec, zapped, plot
  nloop=3                         ; # iterations
  nextra=5d0       ; # pixels expected to be outside sigma factor in image
;This ignores low value pixels near edge that only work to 
; artificially increase sigma and thus make it harder to detect bad
; pixels
  in=where (allspec LT 0 )
  allspec[in]=0
  print, 'negative pixels,', n_elements(in)
  
  yapt1=yapt-4 
  allspec1=allspec[*,4:-5]
 
  inp=1-nextra/4d0/xapt/yapt1
  n=(2d0)^(.5)*inverf(inp) ; sigma rejection factor
  ; Factor chosen such that only 5 pixels in image can be expected
  ; to naturally lie outside the sigma range. I then can correct the pixels
  ; outside the sigma range without worry about overcorrecting
  a=0
  medi=dblarr(5)
  y_medi=findgen(5)*10
  sigma=dblarr(2*xapt)

 
  FOR j = 0, nloop-1 DO BEGIN          ; loop for nloop iterations

; Account for edge effects: If the median of the row is much much
; different then the median of all rows, then don't change the pixels
     row_m=MEDIAN(allspec, Dimension=1)
     row_m=rebin(reform(row_m,1,2*yapt), 2*xapt, 2*yapt)
   ;  med_of_rows=MEDIAN(row_m[4:-5])
   ;  sigma_row=SQRT(TOTAL((row_m[4:-5]-med_of_rows)^2)/(2*yapt1))
   ;  row_med=rebin(reform(row_m,1,2*yapt), 2*xapt, 2*yapt)
   ;  dif_row=ABS(row_med-med_of_rows)
    
; Take care of any other bad pixels by comparing each pixel to the
; median of it's own column
     column_med=MEDIAN(allspec1, Dimension=2)
     column_med1=rebin(column_med, 2*xapt, 2*yapt1)
     sigma=SQRT(TOTAL((allspec1-column_med1)^2,2)/(2*yapt1)) 
    ; s=sigma
     sigma=rebin(sigma, 2*xapt,2*yapt)
     col_med=rebin(column_med, 2*xapt, 2*yapt)
     dif=ABS(allspec-col_med)
    ; column_med=MEDIAN(allspec, Dimension=2)
    ; c2=column_med
    ; column_med=rebin(column_med, 2*xapt, 2*yapt)
    ; sigma=SQRT(TOTAL((allspec-column_med)^2,2)/(2*yapt)) 
    ; STOP, 'check sigma + column_med'
    ; sigma=rebin(sigma, 2*xapt,2*yapt)
    ; dif=ABS(allspec-column_med)
     
; Replace the bad pixels with the median of their column
     xx=indgen(2*yapt)
   ;  rowcut=[max(row_m)*.9]
     rowcut=[median(row_m)*.9]
     rowhigh=[median(row_m)*1.1]
     ycut=rebin(rowcut, 2*xapt, 2*yapt)
     yhigh=rebin(rowhigh, 2*xapt, 2*yapt)
    ; h=poly_fit(xx, row_m, 2)
  ;   p=plot(xx, h[0]+h[1]*xx+h[2]*xx*xx, color='red')
    ;p=plot(xx, row_m[0,*], color='blue',linestyle='', symbol='x')
    ;p=plot(xx, ycut[0, *], color='red', /overplot)
     
     
     index=WHERE(dif GT n*sigma AND row_m GE ycut AND row_m LE yhigh)
     
     ; 2.5 is sufficiently high to exclude edges.Is okay because this 
     ; is a comp. trick, not stat driven
     allspec[index]=col_med[index]
     print, 'bad pixels found =', n_elements(index)
; Plot histogram of pixel distribution for each frame
   ;  IF j EQ 0 AND plot EQ 'on' THEN BEGIN
   ;     FOR x = 0, 2*xapt-1 DO BEGIN      
   ;      hist=HISTOGRAM(allspec[x,*], LOCATIONS=xbin, nbins=150)
   ;      histplot=PLOT(xbin, hist, xrange=[0,max(xbin)],/BUFFER)
   ;      medi[*]=med[x]
   ;      media=PLOT(medi-3*sigma, y_medi, color=1232001, /OVERPLOT)
   ;      media2=PLOT(medi-5*sigma, y_medi, color=112342001, /OVERPLOT)
   ;;      media3=PLOT(medi-6*sigma, y_medi, color=1543132001, /OVERPLOT)
   ;      media4=PLOT(medi,y_medi,color=12433243, /OVERPLOT)
   ;      histplot.save, 'hist.pdf', /APPEND
   ;   ENDFOR
   ;  ENDIF
print,' loop   = ',j+1,'   is finished...',nloop-j,'  to go...'
ENDFOR
 
  zapped=allspec

END

; NAME:
;INVERF -- inverse error function
;     
; PURPOSE:
;       Returns the inverse error function for a given probabality.
;       
;     
; CALLING SEQUENCE:
;       RESULT = INVERF(X)
;     
; INPUTS:
;       X : Number representing the integrated Gaussian probability 
;           distribution between -N standard deviations and +N standard
;           deviations.
;     
; OUTPUTS:
;       Returns the inverse error function of X.  This is the number of
;       standard deviations that the error function needs to be integrated 
;       over in order to yield the value X.
;
; KEYWORDS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS:
;       The error function is normalized to unity, so it is meaningless
;       to ask for the inverse error function of a value larger than
;       one or less than zero.  If you do, the function will stop and 
;       tell you to contemplate your mistake.
;
;       To check behavior visually:
;         IDL> x = dindgen(10001)/10000.
;         IDL> plot, x, abs(x-errorf(inverf(x))), /ylog, yr=[1d-9,1d-1]
;         IDL> oplot, !x.crange, 4.5d-4*[1,1], lines=1
;
;       Hasting's approximation (see NOTES below) is accurate to
;       4.5d-4, and we see this is indeed the case.
;
;       Alternatively, check how large x can become before
;       inverf(errorf(x)) fails, i.e. the difference between
;       x and inverf(errorf(x)) becomes larger than 4.5d-4:
;
;         IDL> x = dindgen(10001)/10000.*6
;         IDL> plot, x, abs(x-inverf(errorf(x))), /ylog, yr=[1d-9,1d-1]
;         IDL> oplot, !x.crange, 4.5d-4*[1,1], lines=1
;
;       So, if x is double precision, inverf(errorf(x)) is accurate
;       up to x=5.34.  If x is floating precision, it is accurate up to
;       x=3.11.  This makes plenty of sense, since beyond this limit out in 
;       the wings of the normal distribution, the error function is exactly 
;       one to machine precision. So the inverse error function can't be
;       expected to differentiate between the error function of 7 and the
;       error function of 200!
;
; NOTES:
;       Using notation from Abramowitz & Stegun, page 931, we see:
;
;         erf(x) = 2*P(x*sqrt(2)) - 1
;         P(x*sqrt(2)) = 1 - Q(x*sqrt(2))
;
;       Given Q(x*sqrt(2))=p, we can use Hasting's approximation 
;       for digital computers (Abramowitz & Stegun, page 933) to get:
;
;         x*sqrt(2) = t - (c0+c1*t+c2*t^2)/(1+d1*t+d2*t^2+d3*t^3) + e(p)
;
;       where |e(p)| < 4.5d-4 and the constants appear in the code below.
;
; MODIFICATION HISTORY:
;       Written Tim Robishaw, Berkeley  22 Feb 2002
;       Stole idea from Carl Heiles.
;       Added machine precision considerations. TR 23 Feb 2002
;-

function inverf, erfx

on_error, 2

; ERROR FUNCTION CAN ONLY RETURN 0 <= ERRORF(X) <= 1, SO
; THE INVERSE ERROR FUNCTION OF ANYTHING GREATER THAN ONE IS
; MEANINGLESS...
if (total((erfx le 1d0) AND (erfx ge 0)) ne N_elements(erfx)) then $
    message, 'The range of the error function is 0 < erf(x) <= 1!'

; IS ERFX WITHIN MACHINE PRECISION OF UNITY...
type = size(erfx,/TYPE)
epsn = (machar(DOUBLE=(type eq 5L))).epsneg

; CHECK THE MACHINE PRECISION HERE...
p = 1d0 - (erfx - epsn*( (erfx - 0.5*epsn) eq fix(1.,TYPE=type) ))
p = 0.5d0*p
t = sqrt(-2d0*alog(p))
num =  (0.010328d0*t + 0.802853d0)*t + 2.515517d0
den = ((0.001308d0*t + 0.189269d0)*t + 1.432788d0)*t + 1d0
return, 0.70710678118654752440d0 * ((t - num/den) > 0)

end; inverf
