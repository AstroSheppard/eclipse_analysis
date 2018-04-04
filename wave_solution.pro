FUNCTION wave_solution, month, dir, PLOTTING=plt, SAVENAME=name, phase=phase, transit=transit


;;;;; Model spectra ;;;;;;;;;;
; Read in model wavelength
  file='./'+month+'/wave.dat'
  OPENR, lun, file, /GET_LUN
  wave=dblarr(8, 152)
  l=8*152
  READF, lun, wave
  wavelength=REFORM(wave,l)
  index=where(wave GT 1000. AND wave LT 1750)
  wavelength=wavelength[index]
  
  FREE_LUN, lun

; Read in model continuum
  file='./'+month+'/continuum.dat'
  OPENR, lun, file, /GET_LUN
  cont=strarr(152)
  READF, lun, cont
  len=10
 

; Convert file to double array
  con=dblarr(l)
  k=0
  FOR i=0, 151 DO BEGIN
     FOR j=0, 7 DO BEGIN
        con[k]=double(STRMID(cont[i], len*j,len))
        k=k+1
     ENDFOR
  ENDFOR
  
  cont=con[index]

  FREE_LUN, lun

; Read in model lines
  file='./'+month+'/kurucz.dat'
  OPENR, lun, file, /GET_LUN
  line=strarr(152)
  READF, lun, line
  len=10
  
  ; Convert to double array
  lines=dblarr(l)
  k=0
  FOR i=0, 151 DO BEGIN
     FOR j=0, 7 DO BEGIN
        lines[k]=double(STRMID(line[i], len*j,len))
        k=k+1
     ENDFOR
  ENDFOR

  line=lines[index]
; Final model
;  ener=6.626e-27*2.998e10/wavelength/1e-8
  f=(line+cont)/wavelength/wavelength


  FREE_LUN, lun

; Read in sensitivity file
  sfile='./wavelength_sol/sens.fits'
  sens=mrdfits(sfile, 1, header)
  wssens=sens.wavelength
  sensitivity=sens.sensitivity
 ; ener=6.626e-27*2.998e10/wssens/1e-7
  through=sensitivity;*ener
; Convert model wavelength to same units (angstroms)
  wavelength=wavelength*10d0
; Interpolate sensitivity array to model wavelength grid
  result=interpol(through, wssens, wavelength)
 
; Convolve model spectrum (f)
 ; inp=dblarr(2, n_elements(index))
 ; inp[0,*]=wavelength
 ; inp[1,*]=f
  ;rpower=[130., 0] ;Resolving power =130, 0 = gaussian kernel
  ;new=convolveinst(inp, rpower)
  ;f=reform(new[1,*])
  ;wavelength=reform(new[0,*])
 


;;;;; Observed Spectrum ;;;;;;;;;;;;


; Find an exposure of just the stellar spectrum (either in eclipse or
; before ingress
  data=FILE_SEARCH('./zapped2017/'+ month +'/'+ dir +'/final/*.zap.fits')
  
  IF ~keyword_set(transit) THEN BEGIN
     IF keyword_set(phase) THEN BEGIN
        IF phase eq 1 then exp=orbits(data, /august1)
        IF phase eq 2 THEN exp=orbits(data, /august2)
     ENDIF ELSE BEGIN
        exp=orbits(data)
     ENDELSE
  ENDIF ELSE BEGIN
     exp=orbits(data, /transit)
  ENDELSE

; Read in data, and normaliza the spectrum
  data=FILE_SEARCH('./zapped2017/'+ month +'/'+ dir +'/final/'+string(exp, format='(I03)')+'.zap.fits')
  spec=MRDFITS(data[0],0, /SILENT)
  spec=total(spec, 2)
  err=sqrt(spec)
  err=err/max(spec)
  spec=spec/max(spec)

 ;;;;; test out removing shape
;  x=findgen(n_elements(spec))
;  coeff=poly_fit(x, spec, 2)
;  y=coeff[2]*x*x + coeff[1]*x + coeff[0]
 ; spec=spec/y
  
  ; Get xcen and ycen for this dataset from the photometry file
  raw=FILE_SEARCH('./'+month+'/*ima.fits')
  header=headfits(raw(0))
  fits_info, raw(0), N_ext=n_ext, /SILENT
  test=MRDFITS(raw[0],1)
  xsize=n_elements(test[*,0])
  ysize=n_elements(test[0,*])

; Isolate the photometric data

  test2=strarr(n_elements(raw))
  FOR i=0, n_elements(raw)-1 DO BEGIN
     exp=MRDFITS(raw[i],0,header,/SILENT)
     test2[i]=fxpar(header,'OBSTYPE')

  ENDFOR
  index=WHERE(test2 NE 'IMAGING       ')
  
  IF index[0] NE -1 THEN BEGIN
     REMOVE, index, raw
     print, n_elements(index), " points have been removed."
  ENDIF
  
  photo='./'+raw[0]
  expo=MRDFITS(photo, 1, header, /SILENT)

; Find the reference pixel
  xcor=fxpar(header, 'LTV1    ')
  ycor=fxpar(header, 'LTV2    ')
  xref=fxpar(header, 'CRPIX1  ')
  yref=fxpar(header, 'CRPIX2  ')
  
  xcen=xref-xcor
  ycen=yref-ycor

; Find range of pixel scales for a given center y
  limits = pixel_scale(ycen)
  scale=mean(limits)
 

;;;;;; FIT MODEL TO DATA ;;;;;;;;
 
  a=[8.95431D3, 9.35925D-2]
  xshift=0.0 
;  norm=1
  p0=[xcen, scale, xshift]

  parinfo = REPLICATE({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0]}, N_ELEMENTS(p0))
  parinfo[*].value=p0
  parinfo[0].fixed=1
  parinfo[1].limited=[1,1]
  parinfo[1].limits=[limits[0], limits[1]]

  fa={flux:spec, err:err, sensitivity:result, star_lines:f, model_wave:wavelength, a:a}


  parameters=MPFIT('solution',functargs=fa,BESTNORM=bestnorm,COVAR=covar,PERROR=perror,PARINFO=parinfo,niter=niter,maxiter=maxiter,status=status,Nfree=nfree)


;;;;;;;;;; Plot results ;;;;;;;;;;;;
  xlen=n_elements(spec)
  data_wave=dblarr(xlen)
  xref=a[0]+a[1]*(0-parameters[0]) 
  FOR i=0, xlen-1 DO BEGIN
     data_wave[i]=(parameters[1]*(i+parameters[2]))+xref
  ENDFOR
  model=f*result
  model=model/max(model)
 
  ;w=wssens/1e10
  ;throughput=sensitivity*
  ;th=throughput/2/max(throughput)

  IF KEYWORD_SET(plt) THEN BEGIN
                                ;  p=plot(data_wave, spec,
                                ;  color='red', title=Month + Dir
                                ;  +'Wavelength Solution')
     p=plot(data_wave/1e4, spec, color='red', xtitle='Wavelength [$\mu$m]', ytitle='Normalized Flux')
     p=plot(wavelength/1e4, model, color='black', /OVERPLOT, xrange=[1,1.8])
     p=plot(wssens/1e4,through/max(through)/2, color='green', /overplot)
     p=plot(wavelength/1e4, f/max(f), color='blue', /OVERPLOT)
     t=text(.28,.24,'Observed Spectrum',  color='red',font_size=11)
     t=text(.28,.20,'Model Spectrum',  color='black', font_size=11)
     t=text(.55,.24,'Model Stellar Flux',  color='blue', font_size=11)
     t=text(.55,.20,'G141 Sensitivity',  color='green', font_size=11)

  
  ENDIF

  IF KEYWORD_SET(name) THEN BEGIN
     SAVE, data_wave, parameters, a, model,wavelength, filename='./wavelength_sol/'+name+'.sav'
     pdffile='./wavelength_sol/'+name+'_wavefit.pdf'
     p.save, pdffile
  ENDIF
  RETURN, data_wave
END

; Function (from WFC3 paper), to determine max and min scale for given ycen
FUNCTION pixel_scale, y
  min=.0028*y+44.68
  max=.0026*y+45.112
  RETURN, [min,max]
END



FUNCTION solution, p, flux=flux, err=err, sensitivity=sens, star_lines=star, model_wave=wave, a=a

; Convert pixels to wavelengths
  xlen=n_elements(flux)
  x=dblarr(xlen)
  xref=a[0]+a[1]*(0-p[0])
  FOR i=0, xlen-1 DO BEGIN
     x[i]=(p[1]*(i+p[2]))+xref
  ENDFOR


; calculate the model (sensitivity function x stellar spectra)
  
  model=star*sens
  model=model/max(model)
  

; Interpolate model to match data wavelength
  theory=interpol(model, wave, x)
 

  RETURN, (flux-theory)/err
  
END

FUNCTION orbits, data, AUGUST1=a1, AUGUST2=a2, transit=transit
  
  orbit=dblarr(1)
  nexposure=n_elements(data)
  date=dblarr(n_elements(data))
  wl=dblarr(n_elements(data))
  FOR img=0, nexposure-1 DO BEGIN
   ; Recover dates
     exposure=MRDFITS(data[img],0,header,/SILENT)
     date[img]=(fxpar(header,'EXPSTART')+fxpar(header,'EXPEND'))/2  
     wl[img]=total(exposure)
  ENDFOR

  FOR i=0, n_elements(date)-2 DO BEGIN
     t=date[i+1]-date[i]
     IF t*86400 GT 1800 THEN orbit=[orbit, i+1] ; 1800s is ~ half an HST orbit
  ENDFOR
  orbit=[orbit,n_elements(date)]
  IF KEYWORD_SET(transit) then exp=orbit[1]-3 ;Proxy for some out of transit exposure.
  IF ~KEYWORD_SET(transit) THEN BEGIN
     nOrbits=n_elements(orbit)-1
     avg_orbit=dblarr(nOrbits)
     dif=dblarr(nOrbits-1)
     i=0
     WHILE i LT nOrbits DO BEGIN
        avg_orbit[i]=median(wl[orbit[i]:orbit[i+1]-1])
        i=i+1
     ENDWHILE
  
     FOR i=0, nOrbits-2 DO BEGIN
        dif[i]=avg_orbit[i+1]-avg_orbit[i]
     ENDFOR

     eclipse=where(dif eq min(dif))+1
     a=sort(dif)
 
     IF KEYWORD_SET(a1) THEN eclipse=a[2]+1
     IF KEYWORD_SET(a2) THEN eclipse=a[1]+1
 ; IF min(dif)/median(wl) LT -.005 THEN eclipse=a[1]+1
     exp=where(wl eq min(wl[orbit[eclipse]:orbit[eclipse+1]-1]))
  ENDIF
  return, exp
  
END

;**************************************************
; convolveinst() function
; Interactive convolution IDL program
; Version 2.0 G. Villanueva - NASA GSFC - Oct/2009
;**************************************************
function convolveinst, $
      p_spec,          $ ; 2D array containing freq. and model data
      respower           ; The resolving power vector (L/dL values), and Kernel and extra values
                         ; 0: Gaussian kernel, 1: triangular, 2: Sinc function, 3: Voigt profile, 4: User-defined kernel

  ; Load the variables
  spec=p_spec
  mfreq=reform(spec(0,*))
  mspec=reform(spec(1,*))
  sz = n_elements(mspec)
  np = n_elements(respower)
  if np gt 1 then kernel = respower(1) else kernel = 0
  if np gt 2 then values = respower(2:np-1) else values = 1d
  
  ;----------------------------------------
  ; Convolution procedure
  ;----------------------------------------
  ; Convolve the spectrum
  fmid = mean(mfreq)
  fdel = abs(median(deriv(mfreq)))
  rpower = respower[0]
  fwhm=1.d*fmid/(rpower*fdel) & if fwhm lt 0.25 then fwhm = 0.25
  fconvs=long(fwhm*10.0)
    
  if fconvs mod 2 ne 0 then fconvs=fconvs+1
  if kernel eq 0 then begin
    ; Gaussian
    fconvx=(dindgen(fconvs) - fconvs/2d)/(fwhm/2.35482d)
    fconvy=exp(-fconvx^2/2.)
  endif else if kernel eq 1 then begin
    ; Triangular
    imid=long(fconvs/2)
    fconvx=(dindgen(fconvs) - imid)
    fconvy=fconvx
    fconvy(0:imid-1)=(fconvx(0:imid-1)+fwhm)/fwhm
    fconvy(imid:fconvs-1)=1.0 - fconvx(imid:fconvs-1)/fwhm
    fconvy(where(fconvy lt 0))=0.0
  endif else if kernel eq 2 then begin
    ; Sinc
    fconvx=(dindgen(fconvs) - fconvs/2d)/fwhm
    fconvy=sin(fconvx*!PI)/(fconvx*!PI)
    fconvy(fconvs/2)=1.0
  endif else if kernel eq 3 then begin
    ; Voigt
    if n_elements(values) eq 0 then wlfac=1d else wlfac=values(0)
    v0 = fmid
    wg = 1d *  (fmid/rpower)                  ; Gauss width [cm-1]
    wl = wlfac*(fmid/2958d) * (0.75d*0.0127d) ; Lorentz width [cm-1]
    fconvs=long(fwhm*100l)
    if fconvs mod 2 eq 0 then fconvs=fconvs+1
    fconvx=dindgen(fconvs)*fdel + v0 - ((fconvs-1d)/2d*fdel)
    fconvy=voigtprof(v0, wg, wl, fconvx)
  endif else if kernel eq 4 then begin
    if n_elements(values) eq 0 then begin
      print,'No kernel given. Stopping ....'
      stop
    endif
    fconvy=values
  endif
    
  ; Convolve the spectrum   
  cspec = convol(mspec, fconvy, total(fconvy), /edge_truncate, center=1)

  ; Return values
  spec_new=spec
  spec_new(1,*)=cspec
  return, spec_new
  
end
