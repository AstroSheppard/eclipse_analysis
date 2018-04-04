; this will read in the june parameters and calculate a weighted
; mean. Then it will output that as a sav file to the directory so
; that aug and apirl can read it in

PRO weight, planet

dir = './WLresults/best_params/constrained/'+planet+'/'
files=FILE_SEARCH(dir+'*.sav')

IF files[-1] eq 'ztime.sav' THEN files=files[0:-2]
IF n_elements(files) eq 2 then begin

   RESTORE, files[0]+'.sav'
   rfinal=final_results ;residuals, depth, err, tcenter, err, inc, err, mpmsr, err
   RESTORE, files[1]+'.sav'
   ffinal=final_results

   tcenter=[rfinal[3],ffinal[3]]
   tcenter_err=[rfinal[4],ffinal[4]]
   
   tres=dblarr(2)

   wmean, tcenter, tcenter_err, tres
ENDIF ELSE BEGIN
   RESTORE, files[0]+'.sav'
   tcenter=final_results[3]
   tcenter_err=final_results[4]
   tres=[tcenter, tcenter_err]
ENDELSE
   

SAVE, filename=dir+'ztime.sav', tres, tcenter, tcenter_err
;print, tres, tres-floor(tres)
END


PRO wmean, val, dval, results
  dw=total(1D/dval^2)
  sum=total(1D*val/dval^2)/dw
  error=sqrt(1D/dw)
  results[0]=sum
  results[1]=error
END
