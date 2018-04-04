; Whitelight curve analysis

.run ./wave_solution_transit.pro
.run ./wave_solution_transit.pro

visit='wasp19/visit01'
direction='reverse'
savename='wasp19b'


;wave=wave_solution(visit, direction , /plotting, savename=savename)

.run ./run.pro
.run ./transit2018.pro
.run ./whitelight_transit2018.pro
.run ./whitelight_transit2018.pro

run, visit, direction, savefile=savename


