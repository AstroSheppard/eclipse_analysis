PRO multi_con, inpfile, planet
; File to find best starting tcenter, inc, and ar*
nVisit=n_elements(inpfile)
FILE_MKDIR, './WLresults/best_params/constrained/'+planet

FOR i=0, nVisit-1 DO BEGIN

   folder = './eclipse_inputs/'
 
   infile=folder+inpfile[i]+'.sav'
   RESTORE, infile
   date=date1
   spectra=allspecextract1
   inputs=props ;[radius planet, tcenter, inc, MpMsR, period, electron/s, depth]
   best_results=results
   t0=double(best_results[3])


  ;initiate arrays
   results=dblarr(9)      ;residuals, depth, depth_err, tcenter, tc_err, inc, 
  ; inc_err, MpMsR (ar), MpMsR_err (ar_err)

   inputs[1]=t0                 ; set tcenter to best value
   inputs[2]=double(inputs[2])
   inputs[3]=double(inputs[3])
  
   
   whitelight_eclipse2017, inputs, date, spectra, first_orbit_size, results
   best_time=results[3]
   time_error=results[4]
   inputs[1]=best_time
   
   starting_params=[best_time, time_error]
      

   whitelight_eclipse2017, inputs, date, spectra, first_orbit_size, results, /NORANDOMT
      
   final_results=results
   
   savfile2='./WLresults/best_params/constrained/'+planet+'/'+inpfile[i]+'.sav'
   SAVE, filename=savfile2, final_results, starting_params
ENDFOR

END

; Program to determine the expected eclipse time
PRO eclipse_time, date, tc
  ; Inputs
  ; date: 1D array of the date of each exposure (MJD)
  ; properties: 1D array containing the last observed eclipse 
  ; and the period. (MJD, days)
  ; NOTE: If transit has been observed but eclipse hasn't, then 
  ; use last observed transit and add period/2
  time=tc[0]
  period=tc[2]
  i=0
  IF date[0] GT time THEN BEGIN
     WHILE (date[0] - time GT period/2d0) DO BEGIN
        time=time+period
        i=i+1
     ENDWHILE
  ENDIF ELSE BEGIN
     WHILE (time - date[0] GT period/2d0) DO BEGIN
        time=time-period
        i=i+1
     ENDWHILE
  ENDELSE
     tc[0]=double(time)
     err=(tc[1]^2-(i*tc[3])^2)^.5
     tc[1]=err
END
