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

    Hamap = fits.open('/home/raul.gonzalez/raul/STARLIGHTv04/'+gal+'_MUSE.lines.fits')[1].data[1]
    signal = fits.open('/home/raul.gonzalez/raul/STARLIGHTv04/SN_maps/s_SII_'+gal+'.fits')[0].data
    
    noise = fits.open('/home/raul.gonzalez/raul/STARLIGHTv04/SN_maps/n_SII_'+gal+'.fits')[0].data

    print('---- Running ADABIN ----')
    s, n, maps = adap_bin(signal, noise, sn_min)
    valores = np.unique(s)
    
    bin_ = 1
    for valor in valores:
        coor = np.where(s == valor)
        s[coor] = bin_
        bin_ += 1

    print('---- '+gal+' segmented in '+str(bin_)+' regions ----')

    path = '/home/raul.gonzalez/raul/STARLIGHTv04/Final_Sample/'+gal+'/'

    os.mkdir('/home/raul.gonzalez/raul/STARLIGHTv04/Final_Sample/'+gal+'/adabin_out_SII/')

    hdur = fits.PrimaryHDU(s)
    hdur.writeto(path+'seg_adabin_SNSII'+str(sn_min)+'_'+gal+'.fits', overwrite=True)

    bin_map = fits.open(path+'seg_adabin_SNSII'+str(sn_min)+'_'+gal+'.fits')[0].data
    cube = fits.open(path+gal+'.resam.fits')[0].data*1e-20
    cube_error = fits.open(path+gal+'.resam.fits')[1].data*1e-20
    cube_hdr = fits.open(path+gal+'.resam.fits')[0].header

    wl = np.array(ascii.read(path+gal+'_150_150.txt')['col1'])
    flag = np.array(ascii.read(path+gal+'_150_150.txt')['col4'])

    f = open(path+gal+'.in')
    f2 = open(path+gal+'_bin'+str(sn_min)+'_merge.in','w')
    line=f.readlines()
    lines=line[0:15]
    lines[0]=str(int(np.max(bin_map)+1))+'\n'
    lines[2] = path+'/adabin_out_SII/'+'\n'
    lines[4] = path+'/adabin_out_SII/'+'\n'
    for l in lines:
        f2.write(l)
    f2.close()
    f.close()

    table = Table(names=[gal+'_bin_1.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m','Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0', gal+'_bin_1'+'.C11.gm.CCM.BN'],dtype=[object, object, object, object, object, object, object, object]) 

    with open(path+'adabin_out_SII/bin_'+gal+'_SN_'+str(sn_min)+'.map', "a") as f:
        f.write ("bin x y\n")
    print('---- Getting regions spectra ----')
    for i in range(1,int(np.max(bin_map)+1)):
        x_ = np.where(bin_map==i)[0]
        y_ = np.where(bin_map==i)[1]
        bin_ = [i]*len(x_)
        if len(x_)>0:
            with open(path+'/adabin_out_SII/bin_'+gal+'_SN_'+str(sn_min)+'.map', "a") as f:
                np.savetxt(f, np.array([bin_, x_, y_]).T, delimiter=' ', fmt="%s")
        
            spec = np.nansum(cube[:,x_,y_],axis=1)
     
            err = np.sqrt(np.nansum(cube_error[:,x_,y_]**2,axis=1))
    
            table.add_row([gal+'_bin_'+(str(i))+'.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m', 'Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0',gal+'_bin_'+(str(i))+'.C11.gm.CCM.BN'])
            np.savetxt(path+'adabin_out_SII/'+gal+'_bin_'+(str(i))+'.txt', np.array([wl, spec, err, flag]).T, delimiter=' ', fmt="%s")

    table.remove_row(0)
    table.write(path+'tabla.txt', format='ascii')
    filenames = [path+gal+'_bin'+str(sn_min)+'_merge.in', path+'tabla.txt']
    with open(path+gal+'_bin'+str(sn_min)+'.in', 'w') as outfile:
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
