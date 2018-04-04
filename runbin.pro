PRO runbin

data=FILE_SEARCH('./eclipse_inputs/*.sav')
; This only gets names, don't worry
n =n_elements(data)
d=strarr(n)
FOR i=0, n-1 DO BEGIN
   l=strlen(data[i])-15
   d[i]=STRMID(data[i],15,l-4)
ENDFOR

FOR i=0, n-1 DO BEGIN
   binsize, d[i], 3
ENDFOR
; d contains all input files for bin
END
