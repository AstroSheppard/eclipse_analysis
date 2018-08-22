PRO wlplot, file

; This procedure reads in the WL results, finds the best fitting
; model, and plots the phase vs systematic-corrected and
; baseline-normalized flux (with errorbars). It overplots the Mandel
; and Agol eclipse model also. 
folder='./paper/wlcurve/'
dir='./WLresults/'
data=file_search(dir+file+'*.sav')
;col=['red', 'red', 'orange', 'orange', 'orange', 'orange', 'green', 'green', 'blue', 'blue']
n=n_elements(data)
;depth=dblarr(10)
nresi=dblarr(n)
ad=dblarr(n)
derr=dblarr(n)
FOR i=0, n-1 DO BEGIN
   RESTORE, data[i]
  
   best=where (w_q eq max(w_q))

   cor=reform(sys_lightcurve[best,*])
   err=reform(sys_lightcurve_err[best,*])
   phase=reform(sys_lightcurve_x[best,*])
   resids=sys_residuals[best,*]*1e6
   
   derr[i]=marg_depth_err*1e6

   model=reform(sys_model[best,*])
   modelx=reform(sys_model_x[best, *])
   xmin=min(phase)-0.05
   xmax=max(phase)+.05


;A-D test
; Make a PDF and CDF of the residuals with minimal binning
   rms=stddev(resids)
   nres=n_elements(resids)
   nresi[i]=nres
  ; pdf1=double(histogram(resids, locations=xbin1, nbins=nres))
  ; pdf1=pdf1/nres
  ; cdf1=total(pdf1, /cumulative)

  ; IF i eq 3 THEN STOP
   ; Put them in order

   res=resids[sort(resids)]
   num=(dblarr(nres)+1.0)/nres
   cdf1=total(num, /cumulative)
   

  
; Get mean, std dev, and peak value of residuals
; case 3
   avg2=mean(res)
   sig2=sqrt(variance(res))

; case 0
   avg=0
   sig=median(err*1e6) ; good guess, but slightly off. Should use 1e6/sqrt(flux), the errors here can be adjusted by the fitting routine

 ; MAKE Continuous one here
; Normalize (see wikipedia)
   dat=(res-avg)/sig
   dat2=(res-avg2)/sig2

   model_dat=findgen(180)/30 - 3
  ; model_dat2=findgen(200)/20 - 5
   model_res = sig*model_dat + avg
   model_res2=sig2*model_dat + avg2
 ;  print, avg2, sig2, sig
; guass_pdf gives a CDF probability for each normalized data point,
; given a gaussian of mean 'avg' and stddev 'sig'
   test=gauss_pdf(dat)
   test2=gauss_pdf(dat2)
   model_test=gauss_pdf(model_dat)
   

; Calculate A-D number
   sum=0
   FOR j=1, nres DO BEGIN
      sum=sum+(2*j-1)*(alog(test[j-1])+alog(1-test[nres-j]))
   ENDFOR
   
   aa=-nres-sum/nres
   ; Adjust for previously unknown mean/sigma 
 ;  aa=aa*(1+4./nres-25./nres/nres)
   ad[i]=aa
   ; save xbin1, cdf1, res, and test for plotting
  

 ;  p=plot(res, cdf1, symbol='o', color='red')
 ;  p=plot(model_res, model_test, /overplot)
 ;  p=plot(model_res2, model_test, color='purple', /overplot)

   SAVE, filename=folder+file_basename(data[i])+'plotstuff.sav', resids, phase, cor, err, model, modelx, cdf1, res, test, model_res, model_test, model_res2
   print, file_basename(data[i])
   ENDFOR

;print, ad
SAVE, filename='./paper/wlcurve/deptherr/'+file+'.sav', derr, ad, nresi, data
SAVE, filename='./paper/wlerror/deptherr/'+file+'.sav', derr, ad, nresi, data

END



  
 
