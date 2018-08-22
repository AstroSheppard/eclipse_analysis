PRO binpaper, visit, bin

dir='./paper/bincurve/'+visit+'/'
f1=dir+'bin'+bin+'.sav'  ; systematic corrected
f2=dir+'norm'+bin+'.sav' ; resids, cdfs, ad number
;f3='aug1r_rawbin'+bin+'.sav' ; raw fluxes

; Read in raw  flux, error (both normalized too) rawflux, rawerr

restore, f2
; residuals, residuals in cdf form, model cdf, eclipse model x, eclipse model flux
; aa, xbin1, cdf1, res, test, resids
restore, f1
;read in systematic reduced light curve, phase, errors
flux=reform(sys_lightcurve[bestfit,*])
phase =reform(sys_lightcurve_x[bestfit,*])
error=reform(sys_lightcurve_err[bestfit,*])
resids_test=reform(sys_residuals[bestfit, *])
modelx=reform(sys_model_x[bestfit, *])
model=reform(sys_model[bestfit,*])

normalizer=min(model)
flux=flux/normalizer
error=error/normalizer
model=model/normalizer
wave=cen/1e4

min=0
max=-1
offset=0.01

p1=errorplot(phase[min:max], rawflux[min:max], rawerr[min:max], color='blue', linestyle='', symbol='o', layout=[1,4,1], ytitle='Normalized Flux', title=' Raw spectral curve, bin ' + string(wave) + 'micron', xrange=[min(phase)-offset,max(phase)+offset],errorbar_color='green')

med_errnorm=median(rawerr)
t=text(.66,.81, 'Representative Error:' , font_size=8)

p2=errorplot(phase[min:max], flux[min:max], error[min:max], color='blue', linestyle='', symbol='o', layout=[1,4,2], /current, xrange=[min(phase)-offset,max(phase)+offset], yrange=[.998, 1.002], ytitle='Normalized Flux', title='Spectral curve with systematics removed', errorbar_color='purple')
p2=plot(modelx, model, /overplot)

med_err=median(error)
t=text(.66,.58, 'Representative Error:', font_size=8)

flat=dblarr(n_elements(resids))
error=error*1e6
r=reform(resids)
p3=errorplot(phase[min:max], r[min:max], error[min:max], color='red', linestyle='', symbol='o', layout=[1,4,3], /current, xtitle='Phase', ytitle='O-M [ppm]', xrange=[min(phase)-offset,max(phase)+offset], errorbar_color='orange', title='Residuals')
p3=plot(phase, flat, /overplot)

mres=mean(resids)
sres=stddev(resids)

t=text(.73,.32, 'Mean of residuals: '+string(mres, format='(D0.2)'), font_size=8)
t=text(.73,.3, 'RMS of residuals: '+string(sres, format='(D0.1)'), font_size=8)

p4=plot(cdf1, res, symbol='o', title='Normality of Residuals', xtitle='CDF', ytitle='Residuals [ppm]', layout=[1,4,4], /current, color='red')
p4=plot(model_test, model_res, /overplot)

t=text(.25,.18, 'RMS/photon: '+string(ratio, format='(D0.2)'), font_size=8)

p4.save, dir+bin+'binpaper.pdf'

p5=plot(cdf1,res,  symbol='o', title='Normality of Residuals', xtitle='CDF', ytitle='Residuals [ppm]', color='red')
p5=plot(model_test, model_res, /overplot)

p5=plot(model_test, model_res2, color='purple', /overplot)

;t=text(.25,.18, 'A-D number: '+string(ratio, format='(D0.2)'), font_size=8)

p5.save, dir + bin + 'normal.pdf'


END
