from astropy.io import fits, ascii
import numpy as np

table = ascii.read('TABLE_MASTER_Final_Sample.tsv')
gals = list(table['galaxy'])
gals=['NGC1365']
path = 'Final_Sample/'
for gal in gals:
    cube = fits.open(path+gal+'/'+gal+'.fits')[1].data
    cube_hdr = fits.open(path+gal+'/'+gal+'.fits')[1].header

    z = float(table['z'][table['galaxy']==gal].value[0])

    wavelength_axis = cube_hdr['NAXIS3']
    step = cube_hdr['CD3_3']
    initial_wavelength = cube_hdr['CRVAL3']
        
    last_wavelength = initial_wavelength + (step * wavelength_axis)
    num_of_elements = int((last_wavelength - initial_wavelength) / step)
    wavelength = np.linspace(initial_wavelength, last_wavelength, num_of_elements, endpoint=True)
    
    

    SII_wl = 6731 * (1 + z)
        
    wl = np.asarray(wavelength)   
    indx_line = (np.abs(wl - SII_wl)).argmin()

    SN=np.zeros((cube.shape[1],cube.shape[2]))
    signal = np.zeros((cube.shape[1], cube.shape[2]))
    noise = np.zeros((cube.shape[1], cube.shape[2]))
    for i in range(0,cube.shape[1]):
        for j in range(0,cube.shape[2]):
            noise[i,j] = np.mean([np.sqrt(np.var(cube[indx_line - 100:indx_line - 35, i, j])), np.sqrt(np.var(cube[indx_line + 35:indx_line + 100, i, j]))])
            signal[i,j] = (np.max(cube[indx_line - 25:indx_line + 25, i, j])-np.mean([np.mean(cube[indx_line - 100:indx_line-35,i,j]), np.mean(cube[indx_line + 35:indx_line + 100,i,j])]))
       
    hdu_new=fits.PrimaryHDU(signal/noise)
    hdu_s = fits.PrimaryHDU(signal)
    hdu_n = fits.PrimaryHDU(noise)
    hdu_new.writeto('sn_maps/'+gal+'_SN_SII_map.fits', overwrite=True)
    hdu_s.writeto('sn_maps/'+gal+'_S_SII_map.fits', overwrite=True)
    hdu_n.writeto('sn_maps/'+gal+'_N_SII_map.fits', overwrite=True)
