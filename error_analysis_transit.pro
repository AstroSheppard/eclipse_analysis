; Before this run @prep (for every visit for the planet). After, run @multiple if there are multiple visits


.run ./wlplot.pro
.run ./wlpaper.pro
.run ./wlerror.pro


file='wasp19a'
planet='wasp19'
nVisit=1
binsize=4
method=3
bin='10'

wlplot, file
wlpaper, file
wlerror, planet

.run ./binsize.pro
.run ./binfit.pro
.run ./bintransit2018.pro
.run ./bin_prep.pro
.run ./binpaper.pro
.run ./get_binerror.pro
.run ./binerror.pro

.run ./binsize.pro
.run ./binfit.pro
.run ./bintransit2018.pro


binsize, file, method, fix=binsize
bin_prep, file
binpaper, file, bin

IF nVisit eq 1 then get_binerror, file, nVisit
IF nVisit eq 1 then binerror, file

