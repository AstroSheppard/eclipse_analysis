PRO bincorrelated, visit, bin

  dir='./paper/bincurve/'+visit+'/'
  f1=dir+'bin'+bin+'.sav'       ; systematic corrected
;  f2=dir+'norm'+bin+'.sav'      ; resids, cdfs, ad number


  restore, f1
; Read in whitelight phase, flux, error (both normalized too)

 
  photons=rawflux_unnorm ;phot=1e6/sqrt(flux)
  resids=reform(sys_residuals[bestfit,*])*1e6
  
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
  
  p=plot(binsize, phot_err, color='red', symbol='o',xtitle='Exposures per bin', ytitle='Error [ppm]', title=visit + ' bin ' + bin+'  correlated noise test')
  p=plot(binsize, rms, color='blue',symbol='x', /overplot)
  t=text(.6,.7, 'Photon Error', color='red')
  t=text(.6,.65, 'RMS', color='blue')
  SAVE, filename='./paper/binerror/correlated/'+visit+bin+'.sav', phot_err, rms, binsize
  p.save, './paper/binerror/correlated/'+visit+bin+'.pdf'

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
