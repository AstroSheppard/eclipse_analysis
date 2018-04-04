PRO bin_prep, file
  dir='./paper/bincurve/'+file
  data=FILE_SEARCH(dir+'*.sav')
 ; IF strcmp(file_basename(data[0]), file) eq 1 THEN data=data[1:-1] 
  num=n_elements(data)

FOR i=0, num-1 DO BEGIN
   restore, data[i]
   title=strmid(file_basename(data[i]), 0, strlen(data[i])-7)

   resids=sys_residuals[bestfit,*]*1e6
   size=(max(resids)-min(resids))/20
   nres=n_elements(resids)
 ;  pdf1=double(histogram(resids, locations=xbin1, nbins=nres))
 ;  pdf1=pdf1/nres
 ;  cdf1=total(pdf1, /cumulative)
      ;;;; AD test ;;;;;
   res=resids[sort(resids)]
   num=(dblarr(nres)+1.0)/nres
   cdf1=total(num,/cumulative)
   

   avg2=mean(res)
   sig2=sqrt(variance(res))
   photon_error=reform(sys_lightcurve_err[bestfit,*])*1e6
  ; STOP
   avg=0
   sig=phot_err
   
   ratio=sig2/sig
;  maxi=max(pdf1)
   dat=(res-avg)/sig
   dat2=(res-avg2)/sig2

   model_dat=findgen(180)/30 - 3
   model_res = sig*model_dat + avg
   model_res2=sig2*model_dat + avg2
   test=gauss_pdf(dat)
   test2=gauss_pdf(dat2)
   model_test=gauss_pdf(model_dat)
   sum=0

  ; If i eq 14 then STOP, avg2, sig2, sig
   FOR j=1, nres DO BEGIN
      sum=sum+(2*j-1)*(alog(test[j-1])+alog(1-test[nres-j]))
   ENDFOR
   aa2=-nres-sum/nres
 ;  aa=aa2*(1+4./nres-25./nres/nres)
   FILE_MKDIR, dir
   name = dir+'/norm'+string(i, format='(i2.2)')+'.sav'
   name2= dir+'/bin'+string(i, format='(i2.2)')+'.sav'
   SAVE, filename=name, aa2, resids, res, test, cdf1  , model_test, test2, model_res, model_res2
   SAVE, filename=name2, sys_lightcurve, sys_lightcurve_x, sys_lightcurve_err, sys_residuals, sys_model, sys_model_x, bestfit, rawflux, rawerr, cen, ratio, rawflux_unnorm, rawerr_unnorm
ENDFOR   
END

;      cd=plot(xbin1, cdf1, symbol='o', title=title+' cdf', xtitle='Residuals [ppm]', ytitle='Dist function', layout=[2,1,1], buffer=1)
 ;     t=plot(res, test2, color='red', /overplot)
 ;     t=text(.12, .81, 'Case 3 = ' + ad2, font_size=10, color='red')
  ;    pd=plot(xbin1, pdf1, histogram=1, xtitle='Residuals [ppm]', ytitle='Number', title=title+' pdf', layout=[2,1,2], /current)
  ;    pdo=plot(res, ytest, /overplot, color='red')
   ;   t=text(.6, .84, 'Mean = ' + av, font_size=10, color='red')
   ;   t=text(.6, .81, 'Sigma = ' + si, font_size=10, color='red')
