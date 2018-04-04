.run ./wlresids.pro ;keyword for transit
.run ./binsize.pro
.run ./binfit.pro
.run ./bintransit2018.pro

visit='wasp19/visit02'
direction='reverse'
file='wasp19a'
method=3


; Get residuals for bin fitting
wlresids, visit, direction, file, /transit
; Run bin size analysis
.run ./binsize.pro
.run ./binfit.pro
.run ./bintransit2018.pro

binsize, file, method, /transit
