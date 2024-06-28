import os
import shutil
from astropy.io import ascii

row = ascii.read('TABLE_MASTER_Final_Sample.tsv')['cubename','galaxy']
cubename, galaxy = list(row['cubename']), list(row['galaxy'])

parent_dir = '/home/hidra2/rgonzalez/Final_Sample/'

for i in range(0,len(cubename)):
    try:
        os.rename(parent_dir+cubename[i]+'.fits', parent_dir+galaxy[i]+'.fits')
    
        path = os.path.join(parent_dir, galaxy[i])
      
        os.mkdir(path)
    
        shutil.move(parent_dir+galaxy[i]+'.fits', parent_dir+galaxy[i]+'/'+galaxy[i]+'.fits')
    except:
