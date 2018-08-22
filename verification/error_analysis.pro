; Before this run @prep (for every visit for the planet). After, run @multiple if there are multiple visits


.run ./wlplot.pro
.run ./wlpaper.pro
.run ./wlerror.pro
.run ./wlcorrelated.pro

file='hatp41'
planet='hatp41'
nVisit=1
binsize=6
method=3
bin='02'

;wlplot, file
;wlpaper, file
;wlerror, planet
;wlcorrelated, file

.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro
.run ./bin_prep.pro
.run ./bincorrelated.pro
.run ./binpaper.pro
.run ./get_binerror.pro
.run ./binerror.pro


.run ./binsize.pro
.run ./binfit.pro
.run ./adjresidual2017.pro


;binsize, file, method, fix=binsize
bin_prep, file
bincorrelated, file, bin
;binpaper, file, bin

;IF nVisit eq 1 then get_binerror, file, nVisit
;IF nVisit eq 1 then binerror, file

