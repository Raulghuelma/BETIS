# HOW TO BIN AN IFS DATACUBE


Sometimes it is necessary to perform adaptive binning to increase the signal-to-noise ratio (S/N) of your data. 
More often than not, when working with IFU data, you may want to enhance the S/N of the nebular emission lines rather than the stellar continuum 
(see BETIS I paper: https://doi.org/10.1051/0004-6361/202348453). Here are the steps to achieve this:

1. Choose your favorite emission line, for instance, [S II]6716. You can apply Equation 1 from the BETIS I paper by running the **make_snr_maps.py** script on your datacube. This will calculate the S/N for each spaxel and return signal, noise, and S/N maps for the [S II]6716 line. The code can be modified to change the emission line of interest and adjust the pseudo-continuum window accordingly.

2. Running **bin_inputs.py** will then bin the datacube as a function of the S/N map obtained in the first step setting a target S/N. For instance, python3 bin_inputs.py -sn 10 -g NGC863 will generate a segmentation map based on the S/N choosen in step 1 with a target of S/N([S II]) = 10, the integrated spectra (and errors) of each bin in a set of .txt files, as well as a .map file where are saved the coordinates of each bin.

3. With the inputs obtained in 2, you can run **STARLIGHT** to perform the simple stellar popularion synthesis to the integrated spectra (see http://www.starlight.ufsc.br/downloads/ to download STARLIGHT and see documentation)

4. Run python3 **bin_inputs.py** -N NGC863 to get the emission lines features. This code will take the **STARLIGHT** outputs and will firstly calculate the nebualar spectra as the difference between the observed integrated spectra and the pupulation models fitted with STARLIGHT, and will save all nebular spectra in a .fits file. Next, using the python **MPFIT** routine, the code will save a .fits file with 7 extensions: EXT0: Flux, EXT1: flux error, EXT2: velocity, EXT3: velocity error, EXT4: gaussian fit parameters (flux, lambda and sigma), EXT5: gaussian fit parameters errors, EXT6: integrated flux. At the same time, each extension contains one slide for every emission line of interest (see the header of the .py). If oyu only want to fit Ha, Hb and equivalent widths, run **fit_bins_EW.py** instead. **THIS STEP CAN TAKE SEVERAL HOURS, DEPENDING OF THE AMOUNT OF INTEGRATED SPECTRA AND LINES OF INTEREST. HIGHLLY RECOMENDABLE TO RUN THE CODE PARALELLY IN A CLUSTER**
