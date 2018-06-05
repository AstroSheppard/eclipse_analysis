*Explanation of each function*

First, data (ima files) must be manually retrieved via FTP (-ip tags, make sure carriage stripping - cr - is off and binary is correct) from Hubble MAST archive and put into a directory named after the planet.
Then, must create 4 files: 
inputs.dat: System parameters from literature
wave.dat, continuum.dat, and kurucz.dat: wavelength, continuum, and line data from kurucz stellar spectra

Now, use programs.
Separate.py: Groups data by visit (but not direction)
%run separate.py [planet]

cr_bkg.py: Isolates data on exposure, sorts by direction, removes cosmic rays and background. Writes data 
            to zapped2017 'planet/visit/direction/cr/'. Requires some manual input (like window size for bkg subtraction).
            Use 'test' as 3rd keyword to see if bkg removal is working as expected (test_zapping)
%run cr_bkg.py [planet] [visit]

fullzap.py: Corrects bad pixels and finds nominal aperture. Writes to zapped2017/.../final
% run fullzap.py [planet] [visit]

Note: Necessary pip installs are given in requirements.txt. Use "pip install -r requirements.txt" 
(preferably in virtualenv) to be able to run. Also, in iPython help(function name) should give more info.
