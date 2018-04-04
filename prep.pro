.run ./wlresids.pro
.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro

visit='kelt1/visit01'
direction='reverse'
file='kelt1r1'
method=3


; Get residuals for bin fitting
wlresids, visit, direction, file
; Run bin size analysis
.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro

binsize, file, method
