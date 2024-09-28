# HOW TO BIN AN IFS DATACUBE


Sometimes it is necessary to perform adaptive binning to increase the signal-to-noise ratio (S/N) of your data. 
More often than not, when working with IFU data, you may want to enhance the S/N of the nebular emission lines rather than the stellar continuum 
(see BETIS I paper: https://doi.org/10.1051/0004-6361/202348453). Here are the steps to achieve this:

1. Choose your favorite emission line, for instance, [S II]6716. You can apply Equation 1 from the BETIS I paper by running the make_snr_maps.py script on your datacube. This will calculate the S/N for each spaxel and return signal, noise, and S/N maps for the [S II]6716 line. The code can be modified to change the emission line of interest and adjust the pseudo-continuum window accordingly.

2. 
