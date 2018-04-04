PRO binfit, inp_file, x, width, center=center, resids=resids, SAVEFILE=sfile, shift=shift

; Note: x = method used: 1, 2, or 3 = Hannah, Avi, or new

;;; Read in exposures and dates to eventually be fit ;;;;

  savfile1='./eclipse_inputs/'+inp_file+'.sav'
  RESTORE, savfile1
  IF x eq 1 THEN BEGIN
     spectra=allspecextract1
     date=date1
  ENDIF
  IF x eq 2 OR x eq 3 THEN BEGIN
     f2='./WLresiduals/'+inp_file+'.sav'

     RESTORE, f2
     spectra=expos
     hst=hstphase
  ; first, cc, date
  ENDIF



 ;;;; Read in best starting depth, best fixed time + inc ;;;;;;;
  savfile3='./WLresults/best_params/'+inp_file+'.sav'
  RESTORE, savfile3
  
  depth_start=final_results[1]
  best_tcenter=final_results[3]
  best_inc=final_results[5]

;;;;;; Read in wavelength solution ;;;;;;;;
  wfile='./wavelength_sol/'+inp_file+'.sav'
 ; wfile='./wavelength_sol/aug1r.sav'
  wavelength=getwave(wfile)

;;;;; Determine bin sized ;;;;;;;;
  index=where(wavelength GE 11000. and wavelength LE 17000.)
  IF ~KEYWORD_SET(shift) THEN shift=0
  size=n_elements(index[shift:-1])
  nbins=(size/width)


;;;;;; Define bins ;;;;;


  start= index[0+shift]
  nexposure=n_elements(spectra[*,0])
  bins=dblarr(nbins,nexposure, width) 
  center=dblarr(nbins)

  FOR i=0, nbins-1 DO BEGIN
  
     n1=start+width*i
     n2=n1+width-1
     center[i]=(wavelength[n1]+wavelength[n2])/2                      
     bin=spectra[*,n1:n2]
     bins[i,*,*]=bin
    
  ENDFOR

 
  inputs=dblarr(7)
  inputs[0]=props[0]
  inputs[1]=best_tcenter
  inputs[2]=best_inc
  inputs[3]=props[3]
  inputs[4]=props[4]
  inputs[5]=1
  inputs[6]=depth_start 
 


  results=dblarr(3)
  r=dblarr(nbins)
  IF ~KEYWORD_SET(resids) THEN resids=dblarr(4)
  depth_array1=dblarr(nbins,2)
  depth_array2=dblarr(nbins,2)
  count=dblarr(nbins)
 

  IF KEYWORD_SET(sfile) THEN BEGIN
     FOR i=0, nbins-1 DO BEGIN
        print, 'bin' + i
        binned_spectra=reform(bins[i,*,*])
        count[i]=total(binned_spectra)
       ; count[i]=median(total(binned_spectra, 2))
        IF x EQ 1 THEN BEGIN
           savename='./spectra_april/models/hannah/' + inp_file  + STRING(i, format='(I02)')
           hannah2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file ,  SAVEFILE=savename
        ENDIF
       
        IF x EQ 2 THEN BEGIN
           savename='./spectra_april/models/resids/' + inp_file + STRING(i, format='(I02)') 
           residual2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file, first, cc, hst , SAVEFILE=savename
        ENDIF
        IF x EQ 3 THEN BEGIN
           cen=center[i]
           savename='./paper/bincurve/' + inp_file + STRING(i, format='(I02)') 
           adjresidual2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file, first, cc, hst , cen, SAVEFILE=savename
        ENDIF
        depth_array2[i,0]=double(results[0])
        depth_array2[i,1]=double(results[1])
        r[i]=double(results[2])
     ENDFOR
  ENDIF ELSE BEGIN
     FOR i=0, nbins-1 DO BEGIN
        binned_spectra=reform(bins[i,*,*])
        count[i]=median(total(binned_spectra, 2))
      
        IF x EQ 1 THEN BEGIN    
           hannah2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file;, /PLOTTING
        ENDIF
        
        IF x EQ 2 THEN BEGIN
           residual2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file, first, cc, hst;, /PLOTTING 
        ENDIF
        IF x EQ 3 THEN BEGIN
           cen=center[i]
           adjresidual2017, inputs, date, binned_spectra, first_orbit_size, results,inp_file, first, cc, hst;, /PLOTTING 
        ENDIF
        depth_array2[i,0]=double(results[0])
        depth_array2[i,1]=double(results[1])
        r[i]=double(results[2])
     ENDFOR
  ENDELSE
  photon_err=1d0/sqrt(count)
  resids[0]=median(r)*1e6
  resids[1]=median(photon_err)*1e6
 ; resid[0,*]=r*1e6
 ; resid[1,*]=photon_err*1e6
     ; Convert to PPM
  error=depth_array2[*,1]*1e6
  depth=depth_array2[*,0]*1e6

;  resids[2]=median(error)
;  resids[3]=median(depth)
  ; replace resids[3] with variance weighted mean
  wdepth=dblarr(2)
  wmean, depth, error, wdepth
  resids[2]=wdepth[1]
  resids[3]=wdepth[0]
  ;resid[2,*]=error
  ;resid[3,*]=depth
  print, 'resids'
  print, r
  print
  print, 'phot'
  print, photon_err
  IF KEYWORD_SET(sfile) THEN BEGIN
     IF x eq 1 THEN BEGIN
        savename2='./spectra_april/'+inp_file + 'hannah'
     ENDIF
     
     IF x eq 2 THEN BEGIN
        savename2='./spectra_april/'+inp_file 
     ENDIF

     IF x eq 3 THEN BEGIN
        savename2='./spectra_april/'+inp_file + 'adj'
     ENDIF

     ratio=r/photon_err
     index=where(ratio LT 1.5)
  ;   center=center[index]
  ;   depth=depth[index]
  ;   error=error[index]
     
     e=errorplot(center[0:-1], depth[0:-1], error[0:-1], linestyle='', symbol='o', color='blue')
    ; save, filename='resids.sav', r, photon_err
     SAVE, filename=savename2+'.sav', wavelength, center, width, nbins, start, depth, error, count
     e.save, savename2+'.pdf'
 
  ENDIF

END

FUNCTION getwave, file
  RESTORE, file
  return, data_wave
END


PRO wmean, val, dval, results
  dw=total(1D/dval^2)
  sum=total(1D*val/dval^2)/dw
  error=sqrt(1D/dw)
  results[0]=sum
  results[1]=error
END
