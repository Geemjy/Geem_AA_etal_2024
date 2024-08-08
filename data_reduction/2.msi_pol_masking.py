#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 
msi_pol_masking.py - Tools for masking the polmask, bad pixel and the backgroud stars (optional).

Written by Geem. J. 2023 - Oct - 13
ksky0422@gmail.com
"""
import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt

from astropy.io import ascii, fits
from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.wcs import WCS
from astroquery.vizier import Vizier
from photutils.centroids import centroid_com
from photutils.aperture import CircularAperture,CircularAnnulus,aperture_photometry
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
import astroscrappy

import warnings
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs.wcs import FITSFixedWarning
warnings.simplefilter('ignore', category=VerifyWarning)
warnings.filterwarnings('ignore', category=UserWarning, append=True)
plt.rcParams['figure.max_open_warning'] = 0
warnings.filterwarnings('ignore', category=FITSFixedWarning) 

########################################
# Where you shoud enter the input      #
########################################

File_path = os.path.join('/Users/judy/Library/CloudStorage/Dropbox/Research/Pirka_MSI/20230929/596') #where the fits file saved
Flat_path = os.path.join('/Users/judy/Library/CloudStorage/Dropbox/Research/Pirka_MSI/20230929/flat') #where the flats are saved
Observatory = 'q33' #Observatory code
Target_name = '596'
tar_type = 'AST' #standard_star = 'STD', asteroids = 'AST'
star_mag = 17 #masking stars brighter than 'star_mag'th (only for tar_type = AST)
star_width = 10 #radius of the masking size of the stars (only for tar_type = AST)
BAD_SIGCLIP = 5 # ??-Sigma clipping for find the bad pixel

#==========================================


# Funtions
def crop(data,row,col,size):
    row_str, row_end = int(row-size), int(row+size)
    col_str, col_end = int(col-size), int(col+size)
    data_cr = data[row_str:row_end,col_str:col_end]
    return data_cr

def skyvalue(data,y0,x0,r_in,r_out,masking=None):
    #To derive the skyvalue 
    if masking is not None:
        masking = masking.astype(bool)
    else:
        masking = np.zeros(np.shape(data))
    # Determine sky and std
    y_in, y_out = int(y0-r_out), int(y0+r_out)
    x_in, x_out = int(x0-r_out), int(x0+r_out)
    #In the case of the targets at the edge of the images
    if y_in < 0:
        y_in = 0
    if y_out > len(data) :
        y_out = len(data)
    if x_in < 0:
        x_in = 0
    if x_out > len(data[0]):
        x_out =  len(data[0])
        
    sky_deriving_area = data[y_in:y_out, x_in:x_out]
    masking = masking[y_in:y_out, x_in:x_out]
    
    new_mask = np.zeros(np.shape(sky_deriving_area))+1
    for yi in range(len(sky_deriving_area)):
        for xi in range(len(sky_deriving_area[0])):
            position = (xi - r_out)**2 + (yi-r_out)**2
            if position < (r_out)**2 and position > r_in**2:
                new_mask[yi, xi] = 0
    new_mask = new_mask.astype(bool)
    mask = new_mask + masking
    
    Sky_region = np.ma.masked_array(sky_deriving_area, mask)
    std = np.ma.std(Sky_region)
    sky = np.ma.median(Sky_region)
    npix = np.shape(sky_deriving_area)[0]*np.shape(sky_deriving_area)[1] - np.sum(mask)
    
    return(sky, std, npix)

def pill_masking(image,x1,x2,y1,y2,Height,target_radi=8):
    x_star_str, x_star_end = x1, x2
    y_star_str, y_star_end = y1, y2
    Masking_image = np.zeros(np.shape(image))
    for yi in range(len(image)):
        for xi in range(len(image)):
            for star in range(len(x_star_end)):
                if star < (len(x_star_end)/2):
                    y_range = (0,116)
                elif star >= (len(x_star_end)/2):
                    y_range = (130,231)
                height = Height[star]
                if x_star_end[star]-x_star_str[star]==0:
                    slope=0
                else:    
                    slope = (y_star_end[star] - y_star_str[star])/(x_star_end[star]-x_star_str[star])
                y_up = slope *xi + y_star_str[star] + height - slope *x_star_str[star]
                y_low = slope *xi + y_star_str[star] - height - slope *x_star_str[star]
                x_str = min(x_star_str[star],x_star_end[star])
                x_end = max(x_star_str[star],x_star_end[star])
                if y_range[0] <= yi <= y_range[1]:
                    if (xi - x_star_str[star])**2 + (yi-y_star_str[star])**2 < (height)**2:
                        Masking_image[yi,xi] = 1
                    if (xi - x_star_end[star])**2 + (yi-y_star_end[star])**2 < (height)**2:
                        Masking_image[yi,xi] = 1    
                    if yi >= y_low-1 and  y_up+1 >= yi and xi > x_str and x_end > xi:
                        Masking_image[yi,xi] = 1      
    return Masking_image   

def delet(header,keyword):
    try: 
        header[keyword]
    except:
        return header
    else:
        del header[keyword]
        return header    
#==========================================

# 2.
#Making the Masking image. 
#The Masking image is masking the background star and the polarization mask area. 
#The position of the backgound stars are queried from astroquery.jplhorizons.Horizons and astroquery.gaia.Gaia
file = sorted(glob(os.path.join(File_path,'msi*.fits')))
if len(file)==0:
    print('There is no msi*.fits file. Check the directory path. Current path:{0}'.format(File_path))
    
magfile = '.csv' 
order = np.arange(0,len(file),4)
_name_print=0
for z in order:
    SET = [file[z],file[z+1], file[z+2], file[z+3]]
    for i in range(0,4):
        RET = SET[i]  
        if os.path.isfile(File_path+'/mask_'+RET.split('/')[-1]):
            #if the masked file exists, continue
            print(File_path+'/mask_'+RET.split('/')[-1] + ' is already exits.')
            continue
        
        #BRING THE IMAGE & ITS HEADER INFO    
        hdul = fits.open(RET)[0]
        header = hdul.header 
        header['NAXIS'] = 2
        header['BZERO'] = 0.
        header['BSCALE'] = 0.
        header['WCSAXIS'] = 2
        image = hdul.data
        if len(np.shape(image))==3:
            image = image[0]
            del header['NAXIS3']
        
        #MAKE THE MASKED IMAGE
        Mask_image = np.zeros(np.shape(image))
        Masking_image_str = np.zeros(np.shape(image))
        if tar_type == 'AST':
            #MASKING THE BACKGROUND STARS
            #Querying the RA, DEC of target based on JD at exposure start
            JD_str = header['MJD-STR'] #JD at exposure start
            obj = Horizons(id=Target_name,id_type = None,location=Observatory,epochs=JD_str)
            eph = obj.ephemerides()
            ra_str,dec_str = eph['RA'][0], eph['DEC'][0]   
            if _name_print == 0:
                print('Target is {0}'.format(eph['targetname'][0]))

            #Querying the RA, DEC of target based on JD at exposure end
            JD_end = header['MJD-END'] #JD at exposure start
            obj = Horizons(id=Target_name,location=Observatory,epochs=JD_str)
            eph = obj.ephemerides()
            ra_end,dec_end = eph['RA'][0], eph['DEC'][0]

            ra_tar = np.mean([ra_str,ra_end])
            dec_tar = np.mean([dec_str,dec_end])

            #Find the background stars's RA,DEC from Gaia
            coord = SkyCoord(ra=ra_tar, dec=dec_tar, unit=(units.degree, units.degree), frame='icrs')
            radi = units.Quantity(0.01, units.deg)
            result = Vizier.query_region(coord,
                            radius = radi,
                            catalog='I/345/gaia2')
            try:
                result[0]
            except IndexError:
                result = []
            else:
                result = result[0]

            RA_star = []
            DEC_star = []
            g_star = []
            if len(result) == 0:
                MASK = Mask_image
                print('No near-bay stars')
            else:
                for i in range(len(result)):
                    if result['Gmag'][i]<=star_mag:
                        RA_star.append(result['RA_ICRS'][i])
                        DEC_star.append(result['DE_ICRS'][i])
                        g_star.append(result['Gmag'][i])

                #Convert the background stars's RA,DEC to pixel coordinate 
                try:
                    data = ascii.read(RET+magfile)
                except FileNotFoundError:
                    print('There is no csv file containing the center of the targets. Run 1.msi_center_init.py first.')
                    exit()
                xo,yo = data['XCENTER'][1]-1,data['YCENTER'][1]-1 #Center of ordinary light
                xe,ye = data['XCENTER'][0]-1,data['YCENTER'][0]-1 # "  of extraordinary light   
                
                crop_size = 20
                torelance = 5
                crop_image = crop(image, yo,xo,crop_size)
                sky_ordi_,sky_std_ordi,area_ordi = skyvalue(crop_image,crop_size,crop_size,
                                                            crop_size/2,crop_size)
                sum_crop_image = crop_image - sky_ordi_
                x1, y1 = centroid_com(sum_crop_image)
                if str(x1) == 'nan':
                    x1,y1 = xo,yo
                g_init = Gaussian2D(amplitude=sum_crop_image[crop_size,crop_size],
                                    x_mean = x1,y_mean=y1,
                                    theta = 0,
                                    bounds={'x_mean':(crop_size-torelance,crop_size+torelance),
                                            'y_mean':(crop_size-torelance,crop_size+torelance),
                                           'x_stddev':(2*gaussian_fwhm_to_sigma,30*gaussian_fwhm_to_sigma),
                                               'y_stddev':(2*gaussian_fwhm_to_sigma,30*gaussian_fwhm_to_sigma)})
                y, x = np.mgrid[:len(crop_image), :len(crop_image[0])]
                fitter = LevMarLSQFitter()
                fitted = fitter(g_init, x, y, sum_crop_image)
                center_x, center_y = fitted.x_mean.value, fitted.y_mean.value
                xo, yo = center_x + (xo-crop_size), center_y + (yo-crop_size) 
                crop_image = crop(image, ye,xe,crop_size)
                sky_extra_,sky_std_extra,area_extra = skyvalue(crop_image,crop_size,crop_size,
                                                              crop_size/2,crop_size)
                sum_crop_image = crop_image - sky_extra_
                x1, y1 = centroid_com(sum_crop_image)
                if str(x1) == 'nan':
                    x1,y1 = xe,ye
                g_init = Gaussian2D(amplitude=sum_crop_image[crop_size,crop_size],
                                    x_mean = x1,y_mean=y1,
                                    theta = 0,
                                    bounds={'x_mean':(crop_size-torelance,crop_size+torelance),
                                            'y_mean':(crop_size-torelance,crop_size+torelance),
                                            'x_stddev':(2*gaussian_fwhm_to_sigma,30*gaussian_fwhm_to_sigma),
                                            'y_stddev':(2*gaussian_fwhm_to_sigma,30*gaussian_fwhm_to_sigma)})
                y, x = np.mgrid[:len(crop_image), :len(crop_image[0])]
                fitter = LevMarLSQFitter()
                fitted = fitter(g_init, x, y, sum_crop_image)
                center_x, center_y = fitted.x_mean.value, fitted.y_mean.value
                xe, ye = center_x + (xe-crop_size), center_y + (ye-crop_size)
                
                #1) FIND THE STARS IN ORDINAY COMPONENT
                header['CRPIX1'], header['CRPIX2'] = xo,yo
                fits.writeto(File_path+'/new_'+RET.split('/')[-1],image,header,overwrite=True)
                ### (X,Y) for the exposure start
                radis_star, X_str, Y_str = [], [], []   
                header['CRVAL1'], header['CRVAL2'] = ra_str,dec_str
                header = delet(header,'CTYPE3')
                header = delet(header,'CRPIX3')
                header = delet(header,'CRVAL3')
                header = delet(header,'CUNIT3')
                header = delet(header,'CDELT3')
                header = delet(header,'CD3_3')
                header = delet(header,'PS3_3')
                header = delet(header,'LTV3')
                header = delet(header,'LTM3_3')
                header = delet(header,'PS3_0')
                header = delet(header,'PV3_1')
                header = delet(header,'PV3_2')
                header = delet(header,'PS3_1')
                header = delet(header,'PS3_2')
                w =WCS(header)
                for i in range(len(RA_star)):
                    x,y = w.wcs_world2pix(RA_star[i],DEC_star[i],0)
                    X_str.append(float(x))
                    Y_str.append(float(y))
                    width = star_width
                    if g_star[i] < 11:
                        radis_star.append(width)
                    elif 11 < g_star[i] < 15:
                        radis_star.append(width*13/15)    
                    elif g_star[i] >= 15:
                        radis_star.append(width* 10/15)
                ### (X,Y) for the exposure end
                X_end,Y_end = [], []         
                header['CRVAL1'], header['CRVAL2'] = ra_end,dec_end     
                w =WCS(header)
                for i in range(len(RA_star)):
                    x,y = w.wcs_world2pix(RA_star[i],DEC_star[i],0)
                    X_end.append(float(x))
                    Y_end.append(float(y))      

                #2) FIND THE STARS IN EXTRAORDINAY COMPONENT
                header['CRPIX1'], header['CRPIX2'] = xe,ye
                hdul.writeto(File_path+'/new_'+RET.split('/')[-1],overwrite=True)
                header['NAXIS'] = 2
                header['EXTEND'] = 'T'
                header['BZERO'] = 0.
                header['BSCALE'] = 0.
                header['WCSAXIS'] = 2
                header = delet(header,'NAXIS3')
                w =WCS(header)


                ### (X,Y) for the exposure start
                header['CRVAL1'], header['CRVAL2'] = ra_str,dec_str    
                header['BZERO'] = 0.
                header['BSCALE'] = 0.
                header['WCSAXIS'] = 2
                w =WCS(header)
                for i in range(0,len(RA_star)):
                    x,y = w.wcs_world2pix(RA_star[i],DEC_star[i],0)
                    X_str.append(float(x))
                    Y_str.append(float(y))        
                    if g_star[i] < 11:
                        radis_star.append(width)
                    elif 11 < g_star[i] < 15:
                        radis_star.append(width*13/15)    
                    elif g_star[i] >= 15:
                        radis_star.append(width*10/15)
                ### (X,Y) for the exposure end    
                header['CRVAL1'], header['CRVAL2'] = ra_end,dec_end     
                w =WCS(header)
                for i in range(len(RA_star)):
                    x,y = w.wcs_world2pix(RA_star[i],DEC_star[i],0)
                    X_end.append(float(x))
                    Y_end.append(float(y))        
                    
                #3)MASKING ALL STARS IN ORDINARY AND EXTRAORDINARY    
                Masking_image_str = pill_masking(image,X_str,X_end,Y_str,Y_end,radis_star)
                plt.imshow(Masking_image_str)
                os.remove(File_path+'/new_'+RET.split('/')[-1])
        
        # MASK THE COSMIC-RAY
        gain = header['GAIN']
        m_LA,cor_image = astroscrappy.detect_cosmics(image,
                                                     gain = gain,
                                                     readnoise = 4.5,
                                                     sigclip=BAD_SIGCLIP)
        tmLA = m_LA.astype(int)
        MASK = Mask_image + Masking_image_str
        MASK[tmLA == 1 ] = 2
        
        
        # MASKING THE POLARIZATION MASK AREA
        # Flat image is used for masking the polarization mask. 
        # Bring the flat whose the retardar angle is same with that of the proccessed image.
        fi = header['FILTER']
        if header['FILTER']=='Rc':
            fil = 'r'
        elif header['FILTER']=='Ic':
            fil ='i'
        elif header['FILTER']=='V':
            fil='v'
        RET_ANG2 = header['RET-ANG2']        
        FLAT = os.path.join(Flat_path,'dome*{0}*_hwp*{1:2.1f}*.fits'.format(fil,RET_ANG2))
        FLAT = glob(FLAT)[0]
        flat = fits.open(FLAT)[0]
        flat0 = flat.data
        if len(np.shape(flat0)) == 3:
            flat0 = flat0[0]
        MASK[flat0 < 0.85] = 1
        MASK[0:17] = 1
        MASK[118:129] = 1
        MASK[231:271] = 1
        MASK[374:385] = 1
        MASK[486:] = 1    
        fits.writeto(File_path+'/mask_'+RET.split('/')[-1],data=MASK,header=header,overwrite=True)
        print(File_path+'/mask_'+RET.split('/')[-1] +' is created.')
        _name_print = _name_print+1
        
