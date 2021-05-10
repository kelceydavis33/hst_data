#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from astropy.coordinates import FK4
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import glob
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
from astropy.coordinates import SkyCoord
from deepCR import deepCR 
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
plt.rcParams['figure.figsize'] = (20, 10)
obj_path = '/data/home/kdavis/hst_data/region.txt'
#obj_path = input('What is the path to the known objects in a .txt file?')
sav_path = '/data/home/kdavis/hst_data/'
#sav_path = input('Directory where the files should be stored:')
file_path = '/data/home/kdavis/data/*.fits'
#file_path= input('What is the path to the files? Include an identifier or asteriks wildcard identifier.')
paths = glob.glob(file_path)
length = len(paths)
thresh = 0.5
#thresh = float(input('What threshold would you like to use? 0.5 reccomended.'))
lost = []
n = 0
bad_im = 0
for path in paths:
    n +=1
    print(f'Processing image {n} of {length}')
    #try:
    image = fits.getdata(path)
#    names = ['raQSO','decQSO']
    dtype = []
    dtype.append( ('raQSO','U20') )
    dtype.append( ('decQSO','U20') )
    gotoqinfo = np.genfromtxt(obj_path,dtype=dtype)
    coords = []
    for i in range(0, gotoqinfo.shape[0]):
        c = SkyCoord(gotoqinfo[i][0], gotoqinfo[i][1], frame = FK4, unit = (u.hourangle, u.deg))
        coords.append(c)
    with fits.open(path) as hdu:
        wcs = WCS(fobj = hdu, header = hdu[1].header)
        pix_x = []
        pix_y = []
        for coord in coords:
            ra = coord.ra
            dec = coord.dec
            y, x = wcs.all_world2pix(ra, dec, 1, adaptive = False, ra_dec_order = True)
            pix_x.append(x)
            pix_y.append(y)
    mdl = deepCR(mask = "ACS-WFC-F606W-2-32", inpaint = "ACS-WFC-F606W-2-32", device = "CPU")
    mask, cleaned_image = mdl.clean(image, threshold = thresh)
    vmax = np.median(cleaned_image)+5*np.std(cleaned_image)
    vmin = np.median(cleaned_image)-5*np.std(cleaned_image)
    dx, dy = 0.005, 0.005

    x = np.arange(np.min(pix_x), np.max(pix_x), dx)
    y = np.arange(np.min(pix_y), np.max(pix_y), dy)
    X, Y = np.meshgrid(x, y)
    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    figs, axs = plt.subplots(1,3, frameon = False)
    axs[0].imshow(cleaned_image, vmin = vmin, vmax = vmax, origin = 'lower', cmap = 'gray', interpolation='bilinear',
                  extent=extent)
    ident = path.split('/')[-1].split('.')[0]
    axs[0].set_title(f'Adjusted')
   # plt.savefig(sav_path + f'{ident}_adjusted.png');
    #plt.figure()
    #vmx = np.median(image)+5*np.std(image)
    #vmn = np.median(image)-5*np.std(image)
    axs[1].imshow(image, vmin = vmin, vmax = vmax, origin = 'lower', cmap = 'gray', interpolation='bilinear',
                                   extent=extent)
    axs[1].set_title(f'Original Image')
    axs[2].imshow(mask, origin = 'lower', cmap = 'gray')
    axs[2].set_title(f'ProbabilityMask', interpolation='bilinear',
                                      extent=extent)
    axs[0].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r', interpolation='bilinear',
                                    extent=extent)
    axs[1].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r', interpolation='bilinear',
                                    extent=extent)
    axs[2].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r', interpolation='bilinear',
                                    extent=extent)
    #axs[3].set_xlim(0,
    #axs[3].set_ylim(0, 
    #axs[3].set_title('Known Stars')
   # plt.show()
    plt.savefig(sav_path + f'{ident}.png')
   # plt.savefig(sav_path + f'{ident}_original');
   # except:
    #    bad_im +=1
     #   lost.append(path.split('/')[-1].split('.')[0])
#print(f'{bad_im} images not processed.The following files were not processed:{lost}')
#file_path = input('What is the absolute path to the data files?')
#thresh = input('What threshold? 0.5 reccomended.')
#image = fits.getdata(file_path)[:512,:512]
#mdl = deepCR(mask="ACS-WFC-F606W-2-32", inpaint="ACS-WFC-F606W-2-32", device="CPU")
#mask, cleaned_image = mdl.clean(image, threshold = thresh)
#prob_mask = mdl.clean(image, binary=False)
#plt.imshow(cleaned_image)
#plt.title('Adjusted Image')
#plt.savefig(f'adjusted_thresh{thresh}.png')
#plt.fig()
#plt.imshow(image)
#plt.title('Original Image')
#plt.savefig('orig.png')

