; Need to make array of constrained (ingress and egress) visit(s), as
; well as an array of all visits 
visits=strarr(1)
planet='wasp19'
;visits=strarr(8)
;convisits=strarr[1]
;convisits[0]='wasp19a'
;convisits[1]='junef'
visits[0]='wasp19a'
;visits[1]='wasp19b'
constrained=1
nCon=0
;nCon=n_elements(convisits)
nVisit=n_elements(visits)


.run ./best_params_single.pro
;.run ./multi_con.pro
;.run ./multi_uncon.pro
.run ./whitelight_eclipse2017.pro
.run ./whitelight_eclipse2017.pro
;.run ./weight.pro

IF nVisit eq 1 and constrained then best_params_single, visits[0], /constrained
IF nVisit eq 1 and ~constrained then best_params_single, visits[0], /uncon

; Using the best constrained visit, save both times
IF nCon gt 0 then multi_con, convisits, planet
; Find weighted mean of times save that
IF nCon gt 0 then weight, planet
; Use that as input to get fixed time for other visits
IF nVisit ne 1 and nCon gt 0 then multi_uncon, visits, planet
IF nVisit ne 1 and nCon gt 0 then multi_uncon, convisits, planet, /best
; If none are constrained, use literature as input for all
IF nVisit ne 1 and nCon eq 0 then multi_uncon, visits, planet, /lit

