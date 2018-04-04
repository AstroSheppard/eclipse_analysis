PRO best_params, inpfile, CONSTRAINED=con, uncon=uncon
; File to find best starting tcenter, inc, and ar*
  folder = './eclipse_inputs/'
 
  infile=folder+inpfile+'.sav'
  RESTORE, infile
  date=date1
  spectra=allspecextract1
  inputs=props ;[radius planet, tcenter, inc, MpMsR, period, electron/s, depth]
  best_results=results
  t0=double(best_results[3])


  ;initiate arrays
  results=dblarr(9) ;residuals, depth, depth_err, tcenter, tc_err, inc, 
  ; inc_err, MpMsR (ar), MpMsR_err (ar_err)

  inputs[1]=t0                  ; set tcenter to best value
  inputs[2]=double(inputs[2])
  inputs[3]=double(inputs[3])
  
  IF KEYWORD_SET(con) THEN BEGIN
     whitelight_eclipse2017, inputs, date, spectra, first_orbit_size, results
     best_time=results[3]
     time_error=results[4]
     inputs[1]=best_time
   
     starting_params=[best_time, time_error]
     
     WLsave='./WLresults/'+inpfile
     whitelight_eclipse2017, inputs, date, spectra, first_orbit_size, results, /NORANDOMT, SAVEFILE=WLsave
                
     final_results=results
     savfile2='./WLresults/best_params/'+inpfile+'.sav'
     SAVE, filename=savfile2, final_results, starting_params
  ENDIF
  IF KEYWORD_SET(uncon) THEN BEGIN
        
     infile2='./WLresults/best_params/'+uncon+'.sav'
     RESTORE, infile2
 

;;;;;;; Read in error on period
     per_file='./eclipse_inputs/period_err/'+inpfile+'.dat'
     OPENR, lun, per_file, /GET_LUN
     per_err=0.0
     READF,lun,per_err
     CLOSE, lun
     FREE_LUN, lun
    

     tc=[tres[0],tres[1],inputs[4],per_err] ;best june time, err, period, err
     eclipse_time, date, tc                 ;comment out for june
     inputs[1]=tc[0]
     tc_err=tc[1]
     
     WLsave='./WLresults/'+inpfile
     whitelight_eclipse2017, inputs, date, spectra, first_orbit_size, results, /FIXTIME, SAVEFILE=WLsave
                            
     final_results=results
     final_results[4]=tc_err
  ;   STOP, results[1], results[3]
     savfile2='./WLresults/best_params/'+inpfile+'.sav'
     SAVE, filename=savfile2, final_results
  ENDIF
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
