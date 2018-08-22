PRO binerror, planet
  dir='./paper/binerror/'
  restore, dir+'scatter/'+planet+'*.sav'
  
  ; by bin
  
 
 ; nvis=n_elements(errbin[0,*])
  nbin=n_elements(errbin[*,0])
  color=3112345*indgen(nbin)
 ; first=1
 ; last=-2
 ; leng=(nbin-first+last+1)*nbin
 ;: a=reform(ratbin[first:last,*], leng, 1)
  FOR i=0, nbin-1 DO BEGIN
     col=color[i]
     rat=reform(ratbin[i,*])
     ad=reform(adbin[i,*])
     IF i eq 0 then p=plot(rat, ad, linestyle='', symbol='x', xtitle='RMS/photon', ytitle='AD number', color=col,layout=[2,1,1], title='Colorized by bin')
     IF i ne 0 then  p=plot(rat, ad, linestyle='', symbol='x', color=col,layout=[2,1,2], /overplot)
 
  ENDFOR


  FOR i=0, nbin-1 DO BEGIN
     col=color[i]
     err=reform(errbin[i,*])
     ad=reform(adbin[i,*])
     IF i eq 0 then p=plot(ad, err, linestyle='', symbol='x', xtitle='AD test', ytitle='Depth error', color=col,layout=[2,1,2], /current, title='Colorized by bin')
     IF i ne 0 then  p=plot(ad, err, linestyle='', symbol='x', color=col,layout=[2,1,2], /overplot)
  ENDFOR
  p.save, dir+planet+'_byBin.pdf'

;SAVE, filename='./scat.sav', ad, rms, err, data, adbin, rmsbin, errbin, databin, rat, ratbin
END
