#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
1.msi_center_init.py - Give the initial point to find the center of the targets

Written by Geem. J. 2023 - Oct - 13
ksky0422@gmail.com
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker
import glob
import os
from astropy.io import fits
import pandas as pd
import sys

########################################
# Where you shoud enter the input      #
########################################

File_path = os.path.join('/Users/judy/Library/CloudStorage/Dropbox/Research/Pirka_MSI/20230929/596') #where the fits file saved

#==========================================


target_direc = os.path.join(File_path)
file = sorted(glob.glob(os.path.join(target_direc,'msi*.fits')))
if len(file)==0:
    print('There is no file. Check if the path is correct. PATH={0}'.format(target_direc))
    sys. exit()

def onclick(event):
    #print ('button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata))
    # Record the x location of the user's click in the global variable and close the figure
    global retval
    retval = (event.xdata, event.ydata)
    plt.close()

for n,fi in enumerate(file):
    df = pd.DataFrame({})
    for t,ray in enumerate(['e','o']):
        fig = plt.figure(figsize=(10,7))
        ax = fig.add_subplot(111)

        hdul = fits.open(fi)[0]
        header = hdul.header
        data = hdul.data
        if np.shape(data)==(1, 512, 512):
            data = data[0]
        # data = data[0:240]
        data_sub = np.copy(data)

        data_sub[17:117] = data_sub[17:117] - np.median(data_sub[17:117])
        data_sub[130:230] = data_sub[130:230] - np.median(data_sub[130:230])
        data_sub[272:373] = data_sub[272:373] - np.median(data_sub[272:373])
        data_sub[386:485] = data_sub[386:485] - np.median(data_sub[386:485])

        ax.imshow(data_sub,vmin=np.percentile(data_sub,10),vmax=np.percentile(data_sub,80),
                  cmap='BuPu')
        saturation = np.copy(data_sub)
        saturation[data<54000] = np.nan
        ax.imshow(saturation,vmin=54000,vmax=55000,cmap='magma')
        # ax.set_ylim(274,486)
        ax.set_ylim(0,236)
        ax.set_xlim(10,len(data[0])-10)

        
        ax.set_title('PRESS {0}-ray!'.format(ray)+fi.split('/')[-1] + ' '*3 + 'HWP={0:2.1f}deg, {1}, {2}'.format(header['RET-ANG2'],
                                                                            header['OBJECT'],
                                                                            header['FILTER'])
                     + ' '*3 
                     +'{0}/{1} ({2}%), {3} sets'.format(n,len(file),int(n/len(file)*100),n//4+1)
                     ,fontsize=12)

        retval = -1

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        # Bring up the figure (and wait)
        plt.show()

        print ('User selected point {0}'.format( retval ))
        XCENTER1 = retval[0]
        YCENTER1 = retval[1]+3

        df = pd.concat([df,
                        pd.DataFrame({'filename':[ray],
                                        'XCENTER':[XCENTER1],
                                        'YCENTER':[YCENTER1]})])
        print(df)
    df = df.round({'XCENTER':2,'YCENTER':2})    
    df.to_csv(fi+'.csv',index=False)
