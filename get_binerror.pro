PRO get_binerror, planet, nVisit

; Program that reads in every light curve, plots it, and makes a
; normality of residuals plot
; Also used to get rms,  depth error, and AD number for each bin
; of each visit

  dir='./paper/bincurve/'
  data=FILE_SEARCH(dir+planet+'*.sav')
  nbin=n_elements(data)/nVisit
  num=n_elements(data)
  rms=dblarr(num)
  ad=dblarr(num)
  err=dblarr(num)
  phot=dblarr(num)
 
FOR i=0, num-1 DO BEGIN
   restore, data[i]
   title=strmid(file_basename(data[i]), 0, strlen(data[i])-7)
   err[i]=marg_error*1e6
   resids=sys_residuals[bestfit,*]*1e6
   rms[i]=stddev(resids)
  
   nres=n_elements(resids)
   pdf1=double(histogram(resids, locations=xbin1, nbins=nres))
   pdf1=pdf1/nres
   cdf1=total(pdf1, /cumulative)
 
      ;;;; Case 3

   res=resids[sort(resids)]
 ;  avg=mean(res)
 ;  sig=sqrt(variance(res))

;; Case 0
   avg=0
   sig=phot_err
  
   maxi=max(pdf1)
   dat=(res-avg)/sig
   inf=where(dat gt 9)
   dat[inf]=8
   test=gauss_pdf(dat)
   sum=0
   FOR j=1, nres DO BEGIN
      sum=sum+(2*j-1)*(alog(test[j-1])+alog(1-test[nres-j]))
   ENDFOR
   
   aa=-nres-sum/nres
 ;  aa=aa*(1+4./nres-25./nres/nres)
   ad[i]=aa
   phot[i]=phot_err
;   ad=strtrim(string(aa),1)

  
 ;     cd=plot(xbin1, cdf1, symbol='o', title=title+' cdf', xtitle='Residuals [ppm]', ytitle='Dist function', buffer=1)
 ;     cd=plot(res, test, color='red', /overplot)
 ;     t=text(.12, .81, 'Case 0 = ' + string(ad[i]), font_size=10, color='red')

;      t=text(.6, .84, 'Mean = 0', font_size=10, color='red')
;;      t=text(.6, .81, 'Sigma = ' + string(sig), font_size=10, color='red')
   ;   print, string(i)
  ;    cd.save, './normality.pdf', /append
   
  print, i

ENDFOR   


adbin=reform(ad, nbin, nVisit)
rmsbin=reform(rms, nbin, nVisit)
errbin=reform(err, nbin, nVisit)
databin=reform(data,nbin,nVisit)
photbin=reform(phot, nbin, nVisit)
rat=rms/phot
ratbin=reform(rat, nbin, nVisit)
;p=plot(rms, err, linestyle='', symbol='x', xtitle='RMS residuals', ytitle='Depth Error')
;p=plot(rms, ad, linestyle='', symbol='x',xtitle='RMS residuals', ytitle='AD stat')
;p=plot(ad, err, linestyle='', symbol='x',xtitle='AD stat', ytitle='Depth Error')


SAVE, filename=dir+'../binerror/scatter/'+planet+'.sav', ad, rms, err, data, adbin, rmsbin, errbin, databin, nbin, nVisit, ratbin, photbin, phot, rat
;cd.save, './normality.pdf', /append, /close

END

