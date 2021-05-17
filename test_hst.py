#!/usr/bin/env python
# coding: utf-8

# In[ ]:
#Importing libraries
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
#Setting figure size, change numbers here to adjust
plt.rcParams['figure.figsize'] = (20, 10)
#Setting the path tolocation of region file that is saved as a text file with no object marker information
#Change the path below to the absolute path to your file or comment it out and un-comment the line below to be
#prompted in an ipython terminal
obj_path = '/data/home/kdavis/hst_data/region.txt'
#obj_path = input('What is the path to the known objects in a .txt file?')
#Setting the path where the .png images should be stored
#Again, replace the path below with your path or change the commented line to be prompted in ipython
sav_path = '/data/home/kdavis/hst_data/'
#sav_path = input('Directory where the files should be stored:')
#Setting the path to the folder containing .fits files you would like to process
#Include a wilcard identifier in your path and again change the commented lines to adjust how you
#enter this information
file_path = '/data/home/kdavis/data/*.fits'
#file_path= input('What is the path to the files? Include an identifier or asteriks wildcard identifier.')
#Using glob to generate a list of filenames
paths = glob.glob(file_path)
#Determine the number of files in the directory
length = len(paths)
#Setting the threshold to use in deepCR. 0.5 is reccomended in the doccumentation but this is the
#parameter you should adjust to change how the images are processed
#Again, there is an option to change this to prompt you for this while running ipython
thresh = 0.5
#thresh = float(input('What threshold would you like to use? 0.5 reccomended.'))
#Creating a parameter to track the number of itterations
n = 0
#Looping over the paths
for path in paths[0:3]:
    #Increasing n by 1
    n +=1
    #printing an indication of the number of images to be processed
    print(f'Processing image {n} of {length}')
    #image is the raw image data from the .fits file
    image = fits.getdata(path)
    #Getting the x and y dimensions of the image
    x_shape = image.shape[0]
    y_shape = image.shape[1]
    #Creating a parameter dtype to store information on the data types of objects in the region file
    dtype = []
    #assuming that these will have RA then DEC coordinates
    dtype.append( ('raQSO','U20') )
    dtype.append( ('decQSO','U20') )
    #Generating data from region file
    gotoqinfo = np.genfromtxt(obj_path,dtype=dtype)
    #Define an empty list to store the SkyCoords from the region file
    coords = []
    #Loop over the objects in the region file
    for i in range(0, gotoqinfo.shape[0]):
        #Create SyCoord for object
        c = SkyCoord(gotoqinfo[i][0], gotoqinfo[i][1], frame = FK4, unit = (u.hourangle, u.deg))
        #Append the SkyCoord to the list
        coords.append(c)
    #Set the Header Data Unit (hdu) to be the hdu from the fits file
    with fits.open(path) as hdu:
        #Create an object wcs to store information about the coordinates from the fits file
        wcs = WCS(fobj = hdu, header = hdu[1].header)
        #Create lists for x an y coordinates of region file objects in this projected coordinates system
        pix_x = []
        pix_y = []
        #Loop over coordinates
        for coord in coords:
            #Grab ra and dec atributes from the coordinates
            ra = coord.ra
            dec = coord.dec
            #Project these ra and dec coordinates to pixel coordinates using wcs
            y, x = wcs.all_world2pix(ra, dec, 1, adaptive = False, ra_dec_order = True)
            #Check to see if these objects are in the frame we plan to plot over
            if (x>0)&(y>0)&(x<x_shape)&(y<y_shape):
                #If they are, add them to our list to plot
                pix_x.append(x)
                pix_y.append(y)
    #Creatinf the model using deepCR. These are all parameters you can change to alter the
    #way that the models are produced. More inputs can be found in the deepCR doccumentation
    mdl = deepCR(mask = "ACS-WFC-F606W-2-32", inpaint = "ACS-WFC-F606W-2-32", device = "CPU")
    #Create a probability mask and a cleaned image using the model above
    #If you are hoping to write out the output of this code to a .fits file, cleaned_image is
    #what you should be writing out.
    mask, cleaned_image = mdl.clean(image, threshold = thresh)
    #Vmax and vmin are used to alter the contrast in the images in matplotlib. This is what you alter when
    #you right click and drag your mouse around in DS9. You can change the contrast by replacing the 5
    #with other numbers
    vmax = np.median(cleaned_image)+5*np.std(cleaned_image)
    vmin = np.median(cleaned_image)-5*np.std(cleaned_image)
    #Extent will be used to define the extent of our imshow plots
    extent = [0, x_shape, 0, y_shape]
    #Defining figs and axs so we can create subplots wihtin one figure  
    figs, axs = plt.subplots(1,3, frameon = False)
    #We put our first plot in the 0 position of the plots we already defined. This will plot our cleaned image.
    axs[0].imshow(cleaned_image, vmin = vmin, vmax = vmax, origin = 'lower', cmap = 'gray', 
                  extent=extent)
    #We create an identifyer to grab the name of the image from the path so it can be added to the filename
    #and plot titles
    ident = path.split('/')[-1].split('.')[0]
    #Adding a title to the plot in the 0 position
    axs[0].set_title(f'Adjusted')
    #PLotting the original image in the first position
    axs[1].imshow(image, vmin = vmin, vmax = vmax, origin = 'lower', cmap = 'gray', 
                                   extent=extent)
    #Adding a title to the original image plot
    axs[1].set_title(f'Original Image')
    #Plotting the probability mask in the 2nd position on the figure
    axs[2].imshow(mask, origin = 'lower', cmap = 'gray', extent = extent)
    #Adding a title to the probability mask 
    axs[2].set_title(f'ProbabilityMask')
    #PLotting a scatterplot of open circles at the position of the known stars on the cleaned image
    axs[0].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r') 
    #PLotting a scatterplot of open circles at the position of the known stars on the original image
    axs[1].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r')
    #PLotting a scatterplot of open circles at the position of the known stars on the probability mask
    axs[2].scatter(pix_x, pix_y, s = 40, facecolors = 'none', edgecolors = 'r')
    #Save the figure to the specified save path, using the image name as an identifyer
    plt.savefig(sav_path + f'{ident}.png')


