PRO binplot

data=file_search('*_binsize.sav')
n=n_elements(data)
name=['WASP-18b Visit 1 Forward', 'WASP-18b Visit 1 Reverse', 'WASP-18b Visit 4 Forward', 'WASP-18b Visit 4 Reverse', 'WASP-18b Visit 5 Forward', 'WASP-18b Visit 5 Reverse', 'WASP-18b Visit 3 Forward', 'WASP-18b Visit 3 Reverse', 'WASP-18b Visit 2 Forward', 'WASP-18b Visit 2 Reverse', 'HAT-P-41b', 'WASP-103b Visit 1 Forward','WASP-103b Visit 1 Reverse','WASP-103b Visit 2 Forward','WASP-103b Visit 2 Reverse','WASP-121b', 'WASP-43b Visit 1 Forward','WASP-43b Visit 1 Reverse','WASP-43b Visit 2 Forward','WASP-43b Visit 2 Reverse','WASP-43b Visit 3 Forward','WASP-43b Visit 3 Reverse','WASP-43b Visit 4 Forward','WASP-43b Visit 4 Reverse','WASP-79b']
For i=0, n-1 DO BEGIN
   RESTORE, data[i]

  p=plot(size, resids, symbol='o', color='red', xtitle='Bin size', /xlog,/ylog,xrange=[3.75,30.75], title='Correlated noise test for ' + name[i]+' Scan', ytitle='RMS')
  p=plot(size, photon, symbol='x', color='blue', /OVERPLOT, /xlog,/ylog)
  t=text(.4,.25,'RMS', color='red')
  t=text(.4,.2,'Photon Noise', color='blue')
  p.save, 'correlate.pdf', /APPEND
ENDFOR
END
