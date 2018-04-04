PRO wlcorrelated, file

  dir='./paper/wlcurve/'
  data=file_search(dir+file+'*')


  restore, data[0]
; Read in whitelight phase, flux, error (both normalized too)
  photons=flux[start:xc] ;phot=1e6/sqrt(flux)
  wlphase=phase
  restore, data[1]
  resids=reform(resids)
  n=n_elements(resids)/2
  phot_err=dblarr(n)
  rms=dblarr(n)
  binsize=findgen(n)+1
  FOR i=1, n DO BEGIN
     photons2=dblarr(n_elements(resids)/i)
     resids2=dblarr(n_elements(resids)/i)
     binflux, photons, photons2, i
     binres, resids, resids2, i
     phot_err[i-1]=1e6/sqrt(median(photons2))
     rms[i-1]=stddev(resids2)
  ENDFOR
  
  p=plot(binsize, phot_err, color='red', symbol='o',xtitle='Exposures per bin', ytitle='Error [ppm]', title=file + ' whitelight correlated noise test')
  p=plot(binsize, rms, color='blue',symbol='x', /overplot)
  t=text(.6,.7, 'Photon Error', color='red')
  t=text(.6,.65, 'RMS', color='blue')
  SAVE, filename='./paper/wlerror/correlated/'+file+'.sav', phot_err, rms, binsize
  p.save, './paper/wlerror/correlated/'+file+'.pdf'

END

PRO binflux, in, out, size
  nbins=n_elements(in)/size
  FOR i=0, nbins-1 DO BEGIN
     start=size*i
     fin=size*(i+1)-1
     out[i]=total(in[start:fin])
  ENDFOR
END
PRO binres, in, out, size
  nbins=n_elements(in)/size
  FOR i=0, nbins-1 DO BEGIN
     start=size*i
     fin=size*(i+1)-1
     out[i]=mean(in[start:fin])
  ENDFOR
END
