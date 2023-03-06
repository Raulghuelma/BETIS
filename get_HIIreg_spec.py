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

def getHII(gal):

    
    path = '/home/raul.gonzalez/raul/STARLIGHTv04/Final_Sample/'+gal+'/'

    os.mkdir('/home/raul.gonzalez/raul/STARLIGHTv04/Final_Sample/'+gal+'/HII_regs/')

    seg_map = fits.open('resultados_'+gal+'_'+gal+'_map_seg.fits')[0].data
    
    cube = fits.open(path+gal+'.resam.fits')[0].data*1e-20
    cube_error = fits.open(path+gal+'.resam.fits')[1].data*1e-20
    cube_hdr = fits.open(path+gal+'.resam.fits')[0].header

    wl = np.array(ascii.read(path+gal+'_150_150.txt')['col1'])
    flag = np.array(ascii.read(path+gal+'_150_150.txt')['col4'])

    f = open(path+gal+'.in')
    f2 = open(path+gal+'_reg_merge.in','w')
    line=f.readlines()
    lines=line[0:15]
    lines[0]=str(int(np.max(seg_map)))+'\n'
    lines[2] = path+'/HII_regs/'+'\n'
    lines[4] = path+'/HII_regs/'+'\n'
    for l in lines:
        f2.write(l)
    f2.close()
    f.close()

    table = Table(names=[gal+'_reg_1.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m','Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0', gal+'_reg_1'+'.C11.gm.CCM.BN'],dtype=[object, object, object, object, object, object, object, object]) 

    with open(path+'HII_regs/reg_'+gal+'.map', "a") as f:
        f.write ("reg x y\n")
    print('---- Getting '+str(int(np.max(seg_map)))+' HII regions spectra ----')
    for i in range(1,int(np.max(seg_map)+1)):
        x_ = np.where(seg_map==i)[0]
        y_ = np.where(seg_map==i)[1]
        bin_ = [i]*len(x_)
        if len(x_)>0:
            with open(path+'/HII_regs/reg_'+gal+'.map', "a") as f:
                np.savetxt(f, np.array([bin_, x_, y_]).T, delimiter=' ', fmt="%s")
        
            spec = np.nansum(cube[:,x_,y_],axis=1)
     
            err = np.sqrt(np.nansum(cube_error[:,x_,y_]**2,axis=1))
    
            table.add_row([gal+'_reg_'+(str(i))+'.txt', 'StCv04.MUSE.config', 'Base.CB07.vall.4m', 'Masks.EmLines.NEW.gm', 'CCM', '0.0', '150.0',gal+'_reg_'+(str(i))+'.C11.gm.CCM.BN'])
            np.savetxt(path+'HII_regs/'+gal+'_reg_'+(str(i))+'.txt', np.array([wl, spec, err, flag]).T, delimiter=' ', fmt="%s")

    table.remove_row(0)
    table.write(path+'tabla.txt', format='ascii')
    filenames = [path+gal+'_reg_merge.in', path+'tabla.txt']
    with open(path+gal+'_reg.in', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
                
          
          
PARSER = argparse.ArgumentParser()
PARSER.add_argument('-g', '--gal', type=str, default=None) 
         
if __name__ == '__main__':
    args = PARSER.parse_args()
    getHII(args.gal)
    elapsed = (time.time() - start_time)
    print(str(timedelta(seconds=elapsed)))
