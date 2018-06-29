PRO binsize, file, method, fix=f

; Procedure to find relationship between bin size and error
  IF ~keyword_set(single) THEN sfile='./paper/bincurve/'+file
  IF keyword_set(f) then begin
     del=file_search(sfile+'*')
     IF n_elements(del) gt 1 or del[0] ne '' THEN FILE_DELETE, del, /recursive
     min=f
     max=f
  ENDIF ELSE BEGIN
     min=4
     max=50
  ENDELSE
  size=dblarr(max-min+1)
  residuals=dblarr(4)
  photon=dblarr(max-min+1)
  resids=dblarr(max-min+1)
  depth_err=dblarr(max-min+1)
  depth=dblarr(max-min+1)
;  sfile='./paper/bincurve/'+file
  FOR i=min, max DO BEGIN
     ;FOR x=1, 3 DO BEGIN
     IF i LT 12 THEN j=i
     IF i GE 12 then j=12+2*(i-12)
     IF keyword_set(f) THEN binfit, file, method, j, center=center, resids=residuals, savefile=sfile
     IF ~keyword_set(f) THEN binfit, file, method, j, center=center, resids=residuals
     resids[i-min]=residuals[0]     ;average rms of residuals per bin
     photon[i-min]=residuals[1]     ;average photon noise per bin
     depth[i-min]=residuals[3]      ; median depth of all bins
     depth_err[i-min]=residuals[2]   ; median depth error for all bins
     size[i-min]=j
  ENDFOR

  res=reform(resids)
  pho=reform(photon)
  ratio=res/pho
  print, ratio

  
  IF ~keyword_set(f) THEN BEGIN

     p=plot(size, resids, symbol='o', color='red', xtitle='Bin size', /xlog,/ylog, title='Correlated Noise test for' + file, ytitle='RMS')
     p=plot(size, photon, symbol='x', color='blue', /OVERPLOT, /xlog,/ylog)
;  t=text(.4,.25,'RMS', color='red')
;  t=text(.4,.2,'Photon Noise', color='blue')
;  print, ''
 ; print, resids
 ; print, photon
  
;  p=plot(center, res, symbol='o', color='blue', sym_filled=1)
;  p=plot(center, pho, symbol='square', color='red', /overplot)
 ; p.save, './paper/binerror/'+file+'.pdf'
;  print, ''
     p.save, file+'_binsize.pdf'

     p1=plot(size, depth, title=file, xtitle='Bin Size', ytitle='Error on Eclipse Depth [PPM]')
     p2=plot(size, depth_err, symbol='x', color='blue', /overplot)
     p2.save, file+'_deptherr.pdf'

     SAVE, filename=file+'_binsize.sav', size, resids, photon, depth, depth_err
  ENDIF

END
