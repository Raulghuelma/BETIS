from adabin.adabin import adap_bin, recon_maps, make_seg_map
from adabin.map_plot import _display_map, recon_plot
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
import warnings
import os
import argparse
warnings.filterwarnings("ignore")
from datetime import timedelta
import time
start_time = time.time()

def run_adabin(sn_min, gal):

    table = ascii.read('TABLE_MASTER_Final_Sample.tsv')
    z = float(table['z'][table['galaxy']==gal].value[0])

    signal = fits.open('/home/rgonzalez/part1/sn_maps/'+gal+'_S_SII_map.fits')[0].data
    
    noise = fits.open('/home/rgonzalez/part1/sn_maps/'+gal+'_N_SII_map.fits')[0].data

    print('---- Running ADABIN ----')
    s, n, maps = adap_bin(signal, noise, sn_min)
    valores = np.unique(s)
    
    bin_ = 1
    for valor in valores:
        coor = np.where(s == valor)
        s[coor] = bin_
        bin_ += 1

    print('---- '+gal+' segmented in '+str(int(bin_))+' regions ----')

    path = '/home/rgonzalez/part1/Final_Sample/'+gal+'/'
    try:
        os.mkdir('/home/rgonzalez/part1/Final_Sample/'+gal+'/adabin_out_SII/')
    except: FileExistsError

    hdur = fits.PrimaryHDU(s)
    hdur.writeto(path+'seg_adabin_SNSII_'+str(int(sn_min))+'_'+gal+'.fits', overwrite=True)

    bin_map = fits.open(path+'seg_adabin_SNSII_'+str(int(sn_min))+'_'+gal+'.fits')[0].data
    cube = fits.open(path+gal+'.resam.fits')[0].data*1e-20
    cube_error = fits.open(path+gal+'.resam.fits')[1].data*1e-20
    cube_hdr = fits.open(path+gal+'.resam.fits')[0].header
    
    wavelength_axis = cube_hdr['NAXIS3']
    step = cube_hdr['CD3_3']
    initial_wavelength = cube_hdr['CRVAL3']
    
    last_wavelength = initial_wavelength + (step * wavelength_axis)
    num_of_elements = int((last_wavelength - initial_wavelength) / step)
    wavelength = np.linspace(initial_wavelength, last_wavelength, num_of_elements, endpoint=True)
    
    restframed_wavelength = wavelength
    number_of_wavelengths = np.floor(max(restframed_wavelength)) - np.ceil(min(restframed_wavelength))
    resampled_wavelength = np.arange(number_of_wavelengths) + np.ceil(restframed_wavelength[0])

    bad = np.array([5453, 5468,
                    5565, 5589,
                    5883, 5903,
                    6294, 6306,
                    4350, 4366,
                    7162, 7333, 
                    9900, 11000]) 

    bad_restframed = bad / (1 + z)
    a = bad_restframed.reshape(7, 2)
    
    flag = np.zeros(len(resampled_wavelength), dtype = np.int8)

    for in_ in range(0, 7):
         flag[(a[in_][0] < resampled_wavelength) & (a[in_][1] > resampled_wavelength)] = 3

    f = open('/home/rgonzalez/part1/Final_Sample/temp.in')
    f2 = open(path+gal+'_bin'+str(int(sn_min))+'_merge.in','w')
    line=f.readlines()
    lines=line[0:15]
    lines[0]=str(int(np.max(bin_map)+1))+'\n'
    lines[2] = path+'/adabin_out_SII/'+'\n'
    lines[4] = path+'/adabin_out_SII/'+'\n'
    for l in lines:
        f2.write(l)
    f2.close()
    f.close()

    table = Table(names=[gal+'_bin_1.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m','Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0', gal+'_bin_1'+'.BN'],dtype=[object, object, object, object, object, object, object, object]) 

    with open(path+'adabin_out_SII/bin_'+gal+'_SN_'+str(int(sn_min))+'.map', "a") as f:
        f.write ("bin x y\n")
    print('---- Getting regions spectra ----')
    for i in range(1,int(np.max(bin_map)+1)):
        x_ = np.where(bin_map==i)[0]
        y_ = np.where(bin_map==i)[1]
        bin_ = [i]*len(x_)
        if len(x_)>0:
            with open(path+'/adabin_out_SII/bin_'+gal+'_SN_'+str(int(sn_min))+'.map', "a") as f:
                np.savetxt(f, np.array([bin_, x_, y_]).T, delimiter=' ', fmt="%s")
        
            spec = np.nansum(cube[:,x_,y_],axis=1)
     
            err = np.sqrt(np.nansum(cube_error[:,x_,y_]**2,axis=1))
    
            table.add_row([gal+'_bin_'+(str(i))+'.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m', 'Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0',gal+'_bin_'+(str(i))+'.BN'])
            np.savetxt(path+'adabin_out_SII/'+gal+'_bin_'+(str(i))+'.txt', np.array([resampled_wavelength, spec, err, flag]).T, delimiter=' ', fmt="%s")

    table.remove_row(0)
    table.write(path+'tabla.txt', format='ascii', overwrite=True)
    filenames = [path+gal+'_bin'+str(int(sn_min))+'_merge.in', path+'tabla.txt']
    with open(path+gal+'_bin'+str(int(sn_min))+'.in', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
                
          
          
PARSER = argparse.ArgumentParser()
PARSER.add_argument('-sn', '--sn_min', type=float, default=None)
PARSER.add_argument('-g', '--gal', type=str, default=None) 
         
if __name__ == '__main__':
    args = PARSER.parse_args()
    run_adabin(args.sn_min, args.gal)
    elapsed = (time.time() - start_time)
    print(str(timedelta(seconds=elapsed)))
