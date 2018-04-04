; Whitelight curve analysis

.run ./wave_solution.pro
.run ./wave_solution.pro

visit='kelt1/visit01'
direction='forward'
savename='kelt1f1'


;wave=wave_solution(visit, direction , /plotting, savename=savename)

.run ./run.pro
.run ./eclipse2017.pro
.run ./whitelight_eclipse2017.pro
.run ./whitelight_eclipse2017.pro

run, visit, direction;, savefile=savename


