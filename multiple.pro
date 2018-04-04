; Final analyses. For multiple visit planets, 

.run ./get_binerror.pro
.run ./binerror.pro

nVisit=2
planet='wasp19'

get_binerror, planet, nVisit
binerror, planet

.run ./best.pro

best, planet
