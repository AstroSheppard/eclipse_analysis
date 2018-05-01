Programs to do whitelight analysis on either eclipses or transits. This is the 2nd group of programs to be run

Need sensitivity.fits from HST website for wavelength calibration. Provided in directory.
Also, occultnl.py and mpfit.py are from other sources. I import them to be used in my functions.

**Order and explanations**
wave_solution.py: Fit observed stellar spectrum to Kurucz stellar model spectrum * grism sensitivity to 
                  determine wavelength of each pixel. Uses mpfit.py. Returns array of physical wavelength for each
                  pixel, and writes to csv's in wave_sol (also saves figures there).
%run wave_solution.py [planet] [visit] [direction] [plotting=True] [transit=False]

whitelight2018.py: Fit lightcurve using a grid of systematic models to determine depth, event time, etc. 
                   Uses mpfit to fit over a grid of 50 models depending on wavelength shift, HST phase, and
                   planetary phase. Returns marginalized depth, center time, and incliantion with errors. 
                   Saves a TON of info to csv's to be handled by multi-indexed pandas dataframes. 
                   
                   1: Anything model dependent (corrected light curves, weights, parameters+errors, etc.) is saved
                   to wl_models_info.csv. Can retrieve by df.loc[visit+direction,desired_index,:] and where transit==transit
                   
                   2: The whitelight data for later plots and measures of success (stddev of resids, rms/photon_err) is
                   saved to wl_data.csv.
This is run by other programs (best_params, preprocess_whitelight.py), and has detailed input info in the file.

preprocess_whitelight.py: Preprocess the observations of 1 visit so we can eventually derive a eclipse/transit spectrum. 
                          Specifically, get best priors, allow user to pick how many orbits surrounding event to use, and cut
                          any suspect data. Save processed_data.csv and system_params.csv, also preprocess_info.csv holds what 
                          orbits and sigma cut to use.
                        
%run program.py [planet] [visit] [direction] [Transit=False]
