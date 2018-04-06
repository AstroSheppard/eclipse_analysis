def data_zapping(month):
	""" Read in all of the fits files, hone in on window around source
	and remove comic rays """
  
	#  dir = '/dosquadequis1/ksheppard/Desktop/2016_research/' 
	###  savename = 'WASP18_WFC3_' 

	# Read in the two columns into two arrays with label and data info
	#READCOL, dir + 'data_files/' + month + '.dat', F='A,D', name, data
	
	###  set_plot, 'x'

	# IMAGE FILES ; Find all image files, store in array "raw". Name of
	# image is stored in "header" array
  	print 'test'
  	raw = FILE_SEARCH('./'+month + '/*ima.fits')
  	header=headfits(raw(0))
  	fits_info, raw(0), N_ext=n_ext, /SILENT
  	test=MRDFITS(raw[0],1)
  	xsize=n_elements(test[*,0])
  	ysize=n_elements(test[0,*])
	# Make sure the data is spectroscopic and SPARS10 (scan, not rapid)
  	test2=strarr(n_elements(raw))
  	test3=strarr(n_elements(raw))
  	quality=strarr(n_elements(raw))
  	FOR i=0, n_elements(raw)-1 DO BEGIN
 		exp=MRDFITS(raw[i],0,header,/SILENT)
     		test2[i]=fxpar(header,'OBSTYPE')
    		test3[i]=fxpar(header,'SAMP_SEQ')
     		quality[i]=fxpar(header,'QUALITY')
  	ENDFOR
  	index=WHERE(test2 NE 'SPECTROSCOPIC ' OR (test3 NE 'SPARS10 ' AND test3 NE 'SPARS25 ') OR quality EQ 'LOCKLOST')
  
  	IF index[0] NE -1 THEN BEGIN
    		REMOVE, index, raw
     	print, n_elements(index), " points have been removed."
  	ENDIF

  	n_forward, n_reverse=0,0
	direction=dblarr(n_elements(raw))

 	FOR i=0, n_elements(raw)-1 DO BEGIN
     		exp=MRDFITS(raw[i],0,header,/SILENT)
     		dire=fxpar(header,'POSTARG2')
     		direction[i]=dire
     		IF dire GT 0 THEN BEGIN
       			n_forward+=1
        		IF n_forward EQ 1 THEN BEGIN
           			forward_img=MRDFITS(raw[i],1)
        		ENDIF	
     		ENDIF ELSE BEGIN
        		n_reverse+=1
     			# Currently unnecessary 
       			IF n_reverse EQ 29 THEN BEGIN ; 
          			reverse_img=MRDFITS(raw[i],1)
       		 	ENDIF
     		ENDELSE
  	ENDFOR
 
  	w1=direction[0]
  	w2=direction[1]
  	rwindow=0
  	fwindow=0
	IF w1 LT 0 THEN rwindow=ceil(2*abs(w1)) else fwindow=ceil(2*abs(w1))
	IF w2 LT 0 THEN rwindow=ceil(2*abs(w2)) else fwindow=ceil(2*abs(w2))
	IF rwindow LE 1 THEN rwindow = fwindow
	IF fwindow LE 1 then fwindow =rwindow
  
  	#STOP,rwindow,  fwindow
	IF n_forward EQ 0 THEN BEGIN
		window, xs=xsize, ys=ysize
		irdisp, reverse_img, low=0
		print, "Please click the bottom-left corner."
		cursor,x,y,/DEVICE
		x1=x
		y1=y

		STOP, "Please type '.cont'"
		Print, "Please click the upper-right corner."
		cursor, x, y,/DEVICE
		x2=x
		y2=y
  	ENDIF ELSE BEGIN
	    	window, xs=xsize, ys=ysize
	    	irdisp, forward_img, low=0
	   	print, "Please click the bottom-left corner."
	   	cursor,x,y,/DEVICE
	   	x1=x
	   	y1=y

	   	STOP, "Please type '.cont'"
	   	Print, "Please click the upper-right corner."
	   	cursor, x, y,/DEVICE
	   	x2=x
	   	y2=y
  	ENDELSE

     	xlen=x2-x1
     	ylen=y2-y1
  	IF n_forward GT 0 THEN BEGIN
     		allimages_f = DBLARR(n_forward, xlen, ylen)
     		allheader_f = STRARR(n_forward, n_elements(header))
  	ENDIF
  	IF n_reverse GT 0 THEN BEGIN
     		allimages_r = DBLARR(n_reverse, xlen, ylen)
     		allheader_r = STRARR(n_reverse, n_elements(header))
	ENDIF


  	f=0
  	r=0
	# Subtract out background
	#;;;;test;;;;
	#bkg, raw[45], rwindow, raw_fit
	#MWRFITS, forward_img*46.7, './test_zapping/raw.fits', /CREATE
	#STOP, 'compare reverse_img[30] with raw_fit[30]'
	#;;;; end test ;;;;

	FOR i = 0, n_elements(raw)-1 DO BEGIN
		# bkg, raw(i), raw_fit
     		print,"image", i
                                #  light[i]=total(mrdfits(raw[i],1))
     		#exposure = MRDFITS(raw[i],0,header,/SILENT)
     		#direction = fxpar(header, 'POSTARG2')
     		IF direction[i] GT 0 THEN BEGIN
      
        		bkg, raw[i], fwindow, raw_fit
          
			img = raw_fit[x1:x2-1,y1:y2-1]
			allimages_f[f,*,*]=img
			header=headfits(raw[i])
			allheader_f[f,*]=header
			f+=1
     		ENDIF ELSE BEGIN

        		bkg, raw[i], rwindow, raw_fit
    
        		img = raw_fit[x1:x2-1,y1:y2-1]
			allimages_r[r,*,*]=img
			header=headfits(raw[i])
			allheader_r[r,*]=header
			r+=1
    		ENDELSE
  
  	ENDFOR

	#	 Remove cosmic rays
  	IF n_reverse GT 0 THEN BEGIN
     		zapped, xlen, ylen, allimages_r, zapped_r, xsize
		rzdir='./zapped2017/'+month + '/reverse/cr/'
		FILE_mkdir, rzdir
		r2dir='./zapped2017/'+month + '/reverse/final/'
		FILE_mkdir, r2dir
 	ENDIF
	IF n_forward GT 0 THEN BEGIN
		zapped, xlen, ylen, allimages_f, zapped_f, xsize
	     	fzdir='./zapped2017/'+month + '/forward/cr/'
	     	FILE_mkdir, fzdir
	     	f2dir='./zapped2017/'+month + '/forward/final/'
	     	FILE_mkdir, f2dir
	ENDIF


	# Write reduced data to a directory
  	FOR k = 0, n_forward-1 DO BEGIN
    		filename = fzdir + STRING(k,format='(I03)') + 'f.zap.fits'
	     	zapp_image = reform(zapped_f(k,*,*)) 
	     	MWRFITS, zapp_image, filename, reform(allheader_f(k,*)), /CREATE 
  	ENDFOR
  
  	FOR k = 0, n_reverse-1 DO BEGIN
     		filename = rzdir + STRING(k,format='(I03)') + 'r.zap.fits'
     		zapp_image = reform(zapped_r(k,*,*)) 
     		MWRFITS, zapp_image, filename, reform(allheader_r(k,*)), /CREATE     
  	ENDFOR
  
  	PRINT, 'Finished zapping ' + month + ' visit'
END

def zapped(xlen, ylen, allspec, allspecnowZap, subarray):

	allspecout=allspec
	nspec=n_elements(allspec(*,0,0))    # number of input spectra
	ny=n_elements(allspec[0,0,*])
	nx=n_elements(allspec[0,*,0])

  	nloop= 0                       # iterations
  	nzap =0                           # sigma rejection factor
 	#size=n_elements(allspec[0,0,*])
  	size = subarray
  	aspec=allspec ;initialize arrays
  	med_rows=allspec
  	row1=dblarr(nspec, ny)
  	FOR j = 0, nloop-1 DO BEGIN          # loop for nloop iterations
     		if j eq 0 then n=8 else n = 6
  
		# Find median of each exposure. Normalize each pixel by that value
		#  med_exp=median(reform(allspec, nspec, npixels), dimension=2)

		#  m=rebin(med_exp, nspec, nx,ny)
		#  aspec=allspec/m
		FOR y=0, ylen-1 DO BEGIN
			rows=median(allspec[*,*,y],dimension=2)

			m=rebin(rows, nspec, nx)
			med_rows[*,*,y]=m
			aspec[*,*,y]=allspec[*,*,y]/m

		ENDFOR
	
		median_row=median(med_rows, dimension=3)

		FOR x = 0, xlen-1 DO BEGIN ;loop over x pixels
		# Find median of each column for each exposure
		# med_column=reform(median(allspec(*,x,*), dimension=3))
		#  m=rebin(med_column, nspec, ny)
		# normalize every pixel by the median of its column
		# aspec=reform(allspec[*,x,*])/m
			FOR y = 0, ylen-1 DO BEGIN ;loop over y pixels
				# Find median of every pixel's time series


				med = MEDIAN(aspec(*,x, y)) 
				sigma = SQRT(TOTAL((aspec(*,x, y)-med)^2)/nspec)   

				FOR z = 0, nspec-1 DO BEGIN ; loop through spectra to clean
					#  xx=indgen(ny)
					#  p=plot(xx, med_rows[z,0,*], color='blue', linestyle='', symbol='x')
					#   if z eq 37 and y eq ymin then stop
					dif = ABS(aspec(z,x,y)-med)

					IF dif GT n*sigma THEN BEGIN ;AND med_rows[z,x,y] GT .95*median_row[z,x] AND med_rows[z,x,y] LT 1.05*median_row[z,x] THEN BEGIN

						allspecout(z,x,y) = med*med_rows[z,x,y]
						# if x GT 30 THEN STOP
						#allspecout[z,x,y]=5e4
						allspec[z,x,y]=med*med_rows[z,x,y]
						nzap=nzap+1
					ENDIF
				ENDFOR
			ENDFOR
		ENDFOR
                               
     		print,' loop   = ',j+1,'   is finished...',nloop-j,'  to go...', size
     		allspec=allspecout
     		print, 'Number Zapped is ', nzap
	ENDFOR
  	allspecnowZap=allspecout
  	RETURN
END

# Input 1 fits file, make array with each extension (3-D: [ext, x, y])
def bkg(raw, size_window, output):
	test=mrdfits(raw,1,/SILENT) 
 	test2=mrdfits(raw,0,header,/SILENT) 
  	corr=fxpar(header, 'UNITCORR')
  	fits_info, raw(0), N_ext=ext, /SILENT
  	xlen=n_elements(test[*,0])
  	ylen=n_elements(test[0,*])
  	nFrames=ext/5
  	frames = dblarr(nFrames, xlen, ylen)
  	times=dblarr(nFrames)
  	frame_diffs = dblarr(nFrames-1, xlen, ylen)
  	count=0
  	for j=1, ext-4, 5 DO BEGIN
    		frames[count,*,*]=mrdfits(raw,j,head,/SILENT)
     		times[count]=fxpar(head,'SAMPTIME')
     		count= count + 1
  	ENDFOR
	#  dq=mrdfits(raw, 3, head,/SILENT)
	#  error=mrdfits(raw, 2, head2,/SILENT)
	#  print, head
	#  print, head2
	#  MWRFITS, dq, './test_zapping/dq.fits',  /CREATE
	#  MWRFITS, error, './test_zapping/error.fits',  /CREATE
	#  STOP


	# Check if the units are in electrons or in electrons/s. If per
	# second, multiply frame by samptime. 
	# Zero out background that is far enough from source to be noise (or
	# is from secondary source)
  	IF corr eq 'OMIT    ' THEN BEGIN
     		# Last frame is for tsamp=0. If brightness not in units of e/s, then
		# this last frame will give large negatives. Set it to 0
     		window=floor(size_window*4)
     		# vary the 8
   		#  window=size_window
     		frames[-1,*,*]=0.0
     		ny=n_elements(frames[0,0,*])
    		# window=40
    		# window=12 ; hatp41
     		for j=0, count-2 DO BEGIN
        		f1=frames[j,*,*]      
        		f2=frames[j+1,*,*]     
        		frame_diffs[j,*,*]=f1-f2
        		mrow=maxrow(frame_diffs[j,*,*]) 
        		IF mrow+window GT ny-1 then window = ny-1-mrow 
        		bg=[[[frame_diffs[j,*,0:mrow-window]]],[[frame_diffs[j,*,mrow+window:-1]]]]
       			med=median(bg)
      			#  print, 'bg', med
        		# Subtract that from whole image
       			frame_diffs[j,*,*]=frame_diffs[j,*,*]-med
        		frame_diffs[j,*,0:mrow-window]=0 
        		frame_diffs[j,*,mrow+window:-1]=0
        		count=count+1
        		file='./test_zapping/d'+string(j)+'.fits'
        		file2='./test_zapping/f'+string(j)+'.fits'
        		MWRFITS, reform(frame_diffs[j,*,*]), file, /CREATE
        		MWRFITS, reform(f1), file2, /CREATE
     		ENDFOR
  	ENDIF
  	IF corr eq 'COMPLETE' THEN BEGIN
    		window=size_window
     		ny=n_elements(frames[0,0,*])
     		for j=0, count-2 DO BEGIN
        		f1=frames[j,*,*]*times[j]
        		f2=frames[j+1,*,*]*times[j+1]
        		frame_diffs[j,*,*]=f1-f2
        		mrow=maxrow(frame_diffs[j,*,*]) 
        		IF mrow+window GT ny-1 then window = ny-1-mrow 
        		# Find median of background noise
        		bg=[[[frame_diffs[j,*,0:mrow-window]]],[[frame_diffs[j,*,mrow+window:-1]]]]
        		med=median(bg)
        		# Subtract that from whole image
        		frame_diffs[j,*,*]=frame_diffs[j,*,*]-med
        		# Zero out outside of mask
        		frame_diffs[j,*,0:mrow-window]=0 ;Change 25 row number for dif data
        		frame_diffs[j,*,mrow+window:-1]=0 ; 23 for april
        		count=count+1
        
      			#  file='./test_zapping/fd'+string(j)+'.fits'
      			#  file2='./test_zapping/ff'+string(j)+'.fits'
      			#  MWRFITS, reform(frame_diffs[j,*,*]), file, /CREATE
      			#  MWRFITS, reform(f1), file2, /CREATE
     			#   file='./test_zapping/b'+string(j)+'.fits'
     			#   MWRFITS, reform(frame_diffs[j,*,*]), file, /CREATE     
     		ENDFOR
	ENDIF
 	# IF ext LT 60 THEN print, raw
  	# Since last WASP79 frame is shifted,treat 2nd to last as last frame
	#mrow=maxrow(ff)
	#  ff[0,*,0:mrow-25]=0
	#  ff[0,*,mrow+25:-1]=0
	#  MWRFITS, reform(ff),'./test_zapping/a6.fits', /CREATE
	#  frame_diffs[-1,*,*]=ff
	#  frame_diffs[-1,*,*]=0
	## Instead of all this, I set the last frame to zero (multiplying it
	## by t=0 effectively does this as well
  
  	output=total(frame_diffs,1)
  
 	# MWRFITS, reform(output), './test_zapping/fbkg.fits', /CREATE

END

def maxrow(frame):
	f=total(frame,2)
  	f1=where(f eq max(f))	
  	return f1
END

     


########################
FUNCTION IRSTR, NUM
;******************************************************************************
;         ** IMPORTANT NOTICE: COPYRIGHTS & ACKNOWLEDGEMENTS ** 
; Copyright (c) 1992, 1993 by: 
; (1) Philip Blanco- CASS, UCSD, 9500 Gilman Drive, La Jolla, CA92093-0111, USA
;     Phone: +1 (619) 534-2943  E-mail: pblanco@ucsd.edu
; (2) Science & Engineering Research Council of the United Kingdom (SERC).
;
; This software may not be copied in whole or in part without the written
; consent of the author. If your research benefits from this code, suitable
; acknowledgement in publications would be appreciated. Bug reports/comments 
; are always welcome. Thanks to R. Pina for ideas and improvements.
;******************************************************************************
;+ IRSTR.PRO      
; Returns a string containing NUM in free-format, with no spaces
; (Note the the GSFC routine STRN does a similar job, and is fancier).
; 
; Parameters (<=input, >=output, !=modified)
; NUM (<) - any number
;
; Calls: none
; 
; History:
; 24 May 92 - written by P. R. Blanco, CASS/UCSD
;
USAGE='<string>=IRSTR(num)'
;-
CASE N_PARAMS() OF 
0: BEGIN
    PRINT, 'Usage: ', USAGE
    RETURN, 0
   END
1: RETURN, STRTRIM(STRING(NUM),2)
ELSE:  BEGIN
        PRINT, 'Incorrect number of arguments- type <var>=IRSTR() for help.'
       RETURN, 0
      END
ENDCASE
;
END

;******************************************************************************
PRO IRDISP, IMAGE,MAG=IMAG,XOFF=DXO,YOFF=DYO,LOW=MINVAL,HIGH=MAXVAL,KEEP=KEEP,quiet=quiet
;******************************************************************************
;         ** IMPORTANT NOTICE: COPYRIGHTS & ACKNOWLEDGEMENTS ** 
; Copyright (c) 1992, 1993 by: 
; (1) Philip Blanco- CASS, UCSD, 9500 Gilman Drive, La Jolla, CA92093-0111, USA
;     Phone: +1 (619) 534-2943  E-mail: pblanco@ucsd.edu
; (2) Science & Engineering Research Council of the United Kingdom (SERC).
;
; This software may not be copied in whole or in part without the written
; consent of the author. If your research benefits from this code, suitable
; acknowledgement in publications would be appreciated. Bug reports/comments 
; are always welcome. Thanks to R. Pina for ideas and improvements.
;******************************************************************************
;+ IRDISP.PRO     
; Produces an image of the "IR_IMAGE" structure, or a 2D array, on the
; current plotting device, magnifying by MAG. 
; The data are plotted between levels MINVAL, MAXVAL. (Irrespective of 
; the values in the data itself, MINVAL and MAXVAL are set
; to the bottom and top colour levels respectively- note this is different
; from TVSCL, which always fills the entire colour table with data). 
; The keyword values are held in a COMMON block IRDISPLAY for use by
; other routines (eg. IRCUR).
;
; Parameters (<=input, >=output, !=modified):
; IMAGE (<) - (2D array, or "IR_IMAGE" structure) - the image to be plotted
; 
; Keywords:
; MAG=IMAG (integer, default 1) - magnification factor for the image display
; XOFF=DXO, YOFF=DYO (integers, default 0) - X and Y offsets (pixels) 
;                    from the origin of the image display.
; LOW=MINVAL (number, default=MIN(image data)) - lowest display level
; HIGH=MAXVAL (number, default=MAX(image data)) - highest display level
; KEEP=keep (logical, default false) - if set, takes the MAX and MIN of
;      the previous display command (IRDISP or IRCONT)
;
; Common blocks
; IRDISP_VALS (>) - holds keyword values for use by IRCUR
;
; Calls: 
; IRSTR ("IR_IMAGE" collection) - returns a number as a string
;
; History:
; 25 May 92  - written by P. R. Blanco, CASS/UCSD
; 24 June 92 - Now the MAG keyword stays in effect for subsequent calls to
;              IRDISP, unless it is changed (i) explicitly in the call by
;              MAG=mag, or (ii) by changing the common block variable IRDMAG
;              at the main level.
; 12 Nov 92  - added keyword KEEP to use value of IRDMIN and IRDMAX currently
;              in IRDISP_VALS common block.
;
COMMON IRDISP_VALS, IRDMAG, IRDXO, IRDYO, IRDMIN, IRDMAX
;
USAGE='IRDISP,<IR_IMAGE image>,[MAG=mag],[XOFF=xo],[YOFF=yo],[LOW=low],' $
     +'[HIGH=high],[/KEEP]'
or_   ='IRDISP, 2D-array, ....'
;-
CASE N_PARAMS() OF
0: BEGIN
    PRINT, 'Usage: ', USAGE
    PRINT, 'or   : ', OR_
    RETURN
   END
1:
;
ELSE: BEGIN
    PRINT, 'Incorrect number of parameters, type IRDISP for help.'
    RETURN
    END
ENDCASE
;
S=SIZE(IMAGE)
NS=N_ELEMENTS(S)
;
; Check to see if IMAGE is a structure, or a 2D array
;
IF (S(NS-2) EQ 8) THEN  BEGIN 
 IMDATA=IMAGE.DATA 
 VALID=WHERE(IMAGE.EXP GT 0.0, N_VAL) 
ENDIF ELSE $
IF (S(0) EQ 2) THEN BEGIN
 IMDATA=IMAGE
 VALID=WHERE(FINITE(IMDATA), N_VAL)
ENDIF ELSE BEGIN
   PRINT, 'Unable to display such an object- wrong dimensions'
   GOTO, ERRSKIP
ENDELSE
;
IF (N_VAL EQ 0) THEN BEGIN
  PRINT, 'No valid pixels in this image'
  goto, errskip
  endif
;
; Get the dimensions and display the scaled image
;
TS=SIZE(IMDATA)
NX=TS(1) & NY=TS(2)
;
; Copy supplied keyword values, or set defaults
;
IF N_ELEMENTS(IMAG) EQ 0 THEN BEGIN
;
; MAG was not specified, use the current value of IRDMAG (if set), or 1
;
   IF N_ELEMENTS(IRDMAG) EQ 0 THEN IRDMAG=1 
ENDIF $
ELSE IRDMAG=FIX(IMAG)
;
IF N_ELEMENTS(DXO) EQ 0 THEN IRDXO=0 ELSE IRDXO=FIX(DXO)
IF N_ELEMENTS(DYO) EQ 0 THEN IRDYO=0 ELSE IRDYO=FIX(DYO)
;
IF NOT KEYWORD_SET(KEEP) THEN BEGIN
;
; Get the max and min values from the data if necessary. Ensure MAX>=MIN
;
 IF N_ELEMENTS(MINVAL) EQ 0 THEN IRDMIN=MIN(IMDATA(VALID)) ELSE IRDMIN=MINVAL
;
 IF N_ELEMENTS(MAXVAL) EQ 0 THEN IRDMAX=MAX(IMDATA(VALID)) ELSE IRDMAX=MAXVAL
;
ENDIF

;
if not  keyword_set(quiet) then PRINT, 'Data has been plotted from LOW='+   $
			   IRSTR(IRDMIN) +' to HIGH='+ IRSTR(IRDMAX)
;
; Convert to byte range for display, using the min amd max limits
; (IRDMIN -->0, and IRDMAX-->!D.N_COLORS in the byte image BTMP).
;
BTMP=BYTSCL(IMDATA,MIN=IRDMIN, MAX=IRDMAX, TOP=!D.N_COLORS)
;
; Display the byte image, suitably magnified and offset from the origin
;
TV, REBIN(BTMP,NX*IRDMAG,NY*IRDMAG, /SAMPLE), IRDXO, IRDYO  
;
RETURN
;
ERRSKIP: RETURN
END
