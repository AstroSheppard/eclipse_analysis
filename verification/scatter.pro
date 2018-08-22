PRO scatter, planet

dat=file_search('./'+planet+'/*.sav')
dir='../../WLresults/'
n=n_elements(dat)
dev=dblarr(n)
rat=dblarr(n)
photerr=dblarr(n)
FOR i=0, n-1 DO BEGIN
   restore, dat[i]
   phot=1e6/median(err)   ;photon error from raw flux
   photerr[i]=phot
   restore, dir+file_basename(dat[i]) ;rms of residuals
   rat[i]=rms/phot
   dev[i]=rms
   STOP
ENDFOR

restore, './deptherr/'+planet+'.sav'
;depth error
;ad test number
;stop
p=plot(rat, ad, linestyle='', symbol='x', xtitle='RMS/photon', ytitle='AD test number', layout=[2,2,1])

p=plot(ad, derr, linestyle='', symbol='x', xtitle='AD test', ytitle='Depth error',layout=[2,2,2], /current)
p=plot(dev, derr, linestyle='', symbol='x', xtitle='RMS', ytitle='Depth error',layout=[2,2,3], /current)
p=plot(rat, derr, linestyle='', symbol='x', xtitle='RMS/photon', ytitle='Depth error',layout=[2,2,4], /current)

p.save, planet+'/wlerror.pdf'


a=sort(rat)
print, dat[a]
print
print, rat[a]
print
print, dev[a]
print
print, photerr[a]
FILE_MKDIR, planet+'/err_info'
save, filename=planet+'/err_info/err_info.sav', rat, dev, photerr, a
END
