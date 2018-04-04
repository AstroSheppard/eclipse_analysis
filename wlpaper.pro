PRO wlpaper, file

dir='./paper/wlcurve/'
data=file_search(dir+file+'*')


restore, data[0]
; Read in whitelight phase, flux, error (both normalized too)
phot=1e6/median(err) 
error=err
wlphase=phase
restore, data[1]
;read in systematic reduced light curve, phase, errors, residuals,
;residuals in cdf form, model cdf, eclipse model x, eclipse model flux
normalizer=min(model)

cor=cor/normalizer
err=err/normalizer
model=model/normalizer

restore, dir+'deptherr/'+file+'.sav'
; Read in A-D number for each curve
index=where(file+'.sav' eq data)
aa=AD[index]


IF file eq 'aug1r' or file eq 'aug1f' then cc=cc-7
p1=errorplot(wlphase[start:cc], fluxnorm[start:cc], errnorm[start:cc], color='blue', linestyle='', symbol='o', layout=[1,4,1], ytitle='Normalized Flux', title=file+' Raw whitelight curve', xrange=[min(phase)-0.03,max(phase)+0.03], errorbar_color='green')

med_errnorm=median(errnorm[start:cc])
t=text(.6,.81, 'Representative Error: ', font_size=8)

p2=errorplot(phase, cor, err, color='blue', linestyle='', symbol='o', layout=[1,4,2], /current, xrange=[min(phase)-0.03,max(phase)+0.03], ytitle='Normalized Flux', title='Whitelight curve with systematics removed', errorbar_color='purple')
p2=plot(modelx, model, /overplot)

med_err=median(err)
t=text(.6,.58, 'Representative Error: ', font_size=8)

flat=dblarr(n_elements(resids))
err=err*1e6
p3=errorplot(phase, reform(resids), err, color='red', linestyle='', symbol='o', layout=[1,4,3], /current, xtitle='Phase', ytitle='Residuals [ppm]', xrange=[min(phase)-0.03,max(phase)+0.03], errorbar_color='orange')
p3=plot(phase, flat, /overplot)

mres=mean(resids)
sres=stddev(resids)

t=text(.6,.32, 'Mean of residuals: '+string(mres, format='(D0.2)'), font_size=8)
t=text(.6,.3, 'RMS of residuals: '+string(sres, format='(D0.1)'), font_size=8)

ratio=sres/phot

p4=plot(cdf1, res, symbol='o', title='Normality of Residuals', xtitle='Dist Function', ytitle='Residuals [ppm]', layout=[1,4,4], /current, color='red', sym_size=1,  yrange=[min(res)-10,max(res)+10])
p4=plot(model_test, model_res, /overplot)
p44=plot(model_test, model_res2, color='purple', /overplot)

t=text(.25,.18, 'RMS/photon = '+string(ratio, format='(D0.2)'), font_size=8)

p4.save, dir+file+'plot.pdf'


p5=plot(cdf1, res, symbol='o', title='Normality of Residuals', xtitle='Dist Function', ytitle='Residuals [ppm]', color='red', sym_size=1.2, yrange=[min(res)-10,max(res)+10])
p5=plot(model_test, model_res, /overplot)

p5=plot(model_test, model_res2, color='purple', /overplot)

;t=text(.25,.18, 'A-D number: '+string(aa, format='(D0.2)'), font_size=8)

p5.save, dir+file+'_normal.pdf'
END
