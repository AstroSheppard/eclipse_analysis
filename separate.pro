PRO separate, dir

; Go to all data file in a directory and separate them based on date
; (so that different visits are separated)

; Read in files
; Find the date of each
; if date2-date1 greater than a hubble orbit period * 2, then save
; files in new directory

  data=FILE_SEARCH('./'+dir+'/'+'*ima.fits')
  nspec=n_elements(data)
  date=dblarr(nspec)
 
; Read in dates for each exposure
  FOR i=0, nspec-1 DO BEGIN
     exp=MRDFITS(data[i], 0, header, /silent)
     date[i]=(fxpar(header,'EXPSTART')+fxpar(header,'EXPEND'))/2.
  ENDFOR

  time=dblarr(nspec-1)
  visit=dblarr(nspec)
  hst_period=95.47
  FOR i=0, nspec-2 DO BEGIN
     t=abs(date[i+1]-date[i])
     t=t*24*60 ; convert to minutes
     time[i]=t
     t=t/hst_period
     IF t GT 3 THEN visit[i]=1 ; If time between data points is greater than 
                                ; 3 HST orbits from previous exposure,
                                ; then I classify it as a new observation
  ENDFOR
 
  nObs=total(visit)+1
  fNames=indgen(nObs)
  numV=0
  dirs='./'+dir+'/visit'+string(fNames,format='(I02)')+'/'
  
  file_mkdir, dirs
  direct=dirs[0]
  FOR i=0, nspec-1 DO BEGIN
   ;  a=mrdfits(data[i],1,/silent)
   ;  header=headfits(data[i])
   ;  filename=direct+string(i, format='(I03)')+'_ima.fits'
   ;  MWRFITS, a, filename, header, /CREATE
     
     FILE_COPY, './'+data[i], direct, /overwrite
     IF visit[i] EQ 1 then begin
        numV=numV+1
        direct=dirs[numV]
     ENDIF
  ENDFOR
END
   
 
