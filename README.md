# HOW TO BIN AN IFS DATACUBE


Sometimes it is necessary to perform adaptive binning to increase the signal-to-noise ratio (S/N) of your data. 
More often than not, when working with IFU data, you may want to enhance the S/N of the nebular emission lines rather than the stellar continuum 
(see BETIS I paper: https://doi.org/10.1051/0004-6361/202348453). Here are the steps to achieve this:

1. Choose your favorite emission line, for instance, [S II]6716. You can apply Equation 1 from the BETIS I paper by running the **make_snr_maps.py** script on your datacube. This will calculate the S/N for each spaxel and return signal, noise, and S/N maps for the [S II]6716 line. The code can be modified to change the emission line of interest and adjust the pseudo-continuum window accordingly.

2. Running **bin_inputs.py** will bin the datacube according to the S/N map obtained in the first step, setting a target S/N. This algorithm makes use of a modified version of the adaptive binning method from Li et al. (2023) (see BETIS I paper) and can be found in the adabin/ directory. For example, python3 bin_inputs.py -sn 10 -g NGC863 will generate a segmentation map based on the S/N chosen in step 1, with a target of S/N([S II]) = 10 (see Figure 2 of the BETIS I paper as an example). This process will also output: the integrated spectra (and errors) for each bin in a set of .txt files, as well as a .map file where are saved the coordinates of each bin.

3. With the inputs obtained in 2, you can run **STARLIGHT** to perform the simple stellar popularion synthesis to the integrated spectra (see http://www.starlight.ufsc.br/downloads/ to download STARLIGHT and see documentation)

4. Run python3 **fit_bin.py** -N NGC863 to get the emission lines features. This code will take the **STARLIGHT** outputs and first calculate the nebular spectra by subtracting the population models fitted with **STARLIGHT** from the observed integrated spectra. It will then save all nebular spectra in a .fits file. Next, using the python **MPFIT** routine (it is important to download the lib/ directory in order to run the codes with all required routines and python packages), the code will save a .fits file with 7 extensions: EXT0: Flux, EXT1: flux error, EXT2: velocity, EXT3: velocity error, EXT4: gaussian fit parameters (flux, lambda and sigma), EXT5: gaussian fit parameters errors, EXT6: integrated flux. At the same time, each extension contains one slide for every emission line of interest (see the header of the .py for details). If you only want to fit H-alpha, H-beta and equivalent  widths, run **fit_bins_EW.py** instead. **THIS STEP CAN TAKE SEVERAL HOURS, DEPENDING OF THE NUMBER OF INTEGRATED SPECTRA AND LINES OF INTEREST. IT IS HIGHLY RECOMENDED TO RUN THE CODE PARALELLY IN A COMPUTING CLUSTER**


5. Finally, with all the fluxes obtained in step 4, you can reconstruct the binned emission line maps by using the coordinates from the .map file or the segmentation map generated in step 2. An example of how to do this is provided in the **recon_maps.ipynb** notebook, where you can follow the steps to recreate the emission line maps based on the segmentation map.

# IMPORTANT CONSIDERATIONS BEFORE BINNING

If you are using the observed and reduced MUSE datacubes (ESO Phase 3), and especially if you plan to perform SSP fitting with **STARLIGHT**, you **MUST** resample the observed spectra of the datacube to match the SSP spectral sampling and apply a correction for Milky Way extinction.

To do this, run the following command:

python3 **resamp_cube.py** -z (redshift) -a (a_v) -n (survey name e.g. MUSE, MANGA or PMAS) -g (Galaxy name)
This will generate a galaxy.resam.fits file, which contains the observed datacube that has been rest-framed, resampled to 1 Å, and corrected for Milky Way extinction.

If you don't want to bin the datacube, you can run **prestarcube.py** instead of **resamp_cube.py** (with the same input commands) to perform the same resampling **AND ALSO** prepare the inputs for running **STARLIGHT** spaxel-by-spaxel. Then, running **fit_emi_abs_ccorr_IFU.py** -N galaxy_name will generate a nebular emission line cube, as well as the same outputs as in step 4, but in the form of emission line maps calculated spaxel-by-spaxel.

![seg_map](https://github.com/user-attachments/assets/5fc49c3b-3bc8-4610-97fe-83de852678f9)


