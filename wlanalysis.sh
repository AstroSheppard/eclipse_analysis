; Whitelight curve analysis

visit='kelt1/visit01'
direction='forward'
savename='kelt1f1'

; Run wave_solution.py

wave=wave_solution(visit, direction , /plotting, savename=savename)
;Run run.py, which should import eclipse2017.py, whitelight_eclipse2017.py

run, visit, direction, savefile=savename


