.run ./wlresids.pro
.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro

visit='wasp19/visit02'
direction='reverse'
file='wasp19a'
method=3


; Get residuals for bin fitting
wlresids, visit, direction, file
; Run bin size analysis
.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro

;binsize, file, method
