def best(planet):

dir='./spectra_april/'
data=file_search(dir+planet+'*.sav')
restore, data[0]
nBin=n_elements(depth)
nVisit=n_elements(data)
margdepths=dblarr(nVisit)
margerrs=dblarr(nVisit)
spectra=dblarr(nVisit, nBin, 2)
basedir='./WLresults/'
FOR i=0, nVisit-1 DO BEGIN

   restore, data[i]
   spectra[i,*,0]=depth
   spectra[i,*,1]=error
   IF i eq 0 then wcenter=center
   WLfile=basedir+file_basename(data[i])
   restore, WLfile
   margdepths[i]=marg_depth
   margerrs[i]=marg_depth_err
ENDFOR
   
  
margWL=dblarr(2)

wmean, margdepths, margerrs, margWL
dif=(-1.0*margdepths+margWL[0])*1e6

FOR i=0, nVisit-1 DO BEGIN
   spectra[i,*,0]=spectra[i,*,0]+dif[i]
ENDFOR

final_spectra=dblarr(nBin,2)
FOR i=0, nBin-1 DO BEGIN
   bindepth=reform(spectra[*,i,0])
   binerror=reform(spectra[*,i,1])
   final_bin=dblarr(2)
   wmean, bindepth, binerror, final_bin
   final_spectra[i,*]=final_bin
ENDFOR

center=wcenter
depth=final_spectra[*,0]
error=final_spectra[*,1]

first=0
last=-1
center=center[first:last]
depth=depth[first:last]
error=error[first:last]
e=errorplot(wcenter, depth, error, linestyle='', symbol='o',color='red',errorbar_color='red', title='Emission Spectra')
e.save, dir+planet+'.pdf'

SAVE, filename=dir+planet+'.sav', center, results, depth, error, dif

;e.close
END


PRO wmean, val, dval, results
  dw=total(1D/dval^2)
  sum=total(1D*val/dval^2)/dw
  error=sqrt(1D/dw)
  results[0]=sum
  results[1]=error
END


