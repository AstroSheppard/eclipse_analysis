*Explanation of each function*

First, data (ima files) must be manually retrieved via FTP (-ip tags, make sure carriage stripping - cr - is off and binary is correct) from Hubble MAST archive and put into a directory named after the planet.
Then, must create 4 files: 
inputs.dat: System parameters from literature
wave.dat, continuum.dat, and kurucz.dat: wavelength, continuum, and line data from kurucz stellar spectra

Now, use programs.
Separate.py: Groups data by visit (but not direction)
%run separate.py [planet]

bkg.py: Isolates data on exposure, sorts by direction, gets initial aperture, and removes background. Writes data 
            to zapped2017 'planet/visit/direction/bkg/'. Requires some manual input (like window size for bkg subtraction).
            Use 'test' as 3rd keyword to see if bkg removal is working as expected (test_zapping)
%run bkg.py [planet] [visit] [test = false]

reduction.py: Finds wavelength solution, removes flat field, masks bad pixels from dq array, removes cosmic rays, sets final aperture, 
              updates errors for background noise, and writes full fits file (reduced data, errors, raw image) to zapped2017/.../final
% run reduction.py [planet] [visit] [direction] [transit=0]

fullzap.py: Contains functions used in reduction.py

sensitivity.fits and flat.fits: From instrument, used in wavelength solution and flat field removal, respectively.

Note: Necessary pip installs are given in requirements.txt. Use "pip install -r requirements.txt" 
(preferably in virtualenv) to be able to run. Also, in iPython help(function name) should give more info (from docstrings).
