PRO run, month, direction, SAVEFILE=name



 ; set_plot, 'x'


; inputs: Epoch time, first orbit, last orbit, sigma, nPasses
  inputs=dblarr(5)
  results=dblarr(9)

  x=0
  y=0
  ;inputs=[2,7,5,1,]
; Run once to get user inputs
  eclipse2017, month, direction, X=x, Y=y, best_results,inputs, /PLOTON, /CHECK,savefile=name, /SAVEDATA


END




