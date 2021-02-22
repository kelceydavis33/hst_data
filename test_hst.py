#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt
import glob
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
from astropy import units as u
from astropy.coordinates import SkyCoord
from deepCR import deepCR as dcr
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
file_path = input('What is the absolute path to the data files?')
thresh = input('What threshold? 0.5 reccomended.')
image = fits.getdata(file_path)[:512,:512]
mdl = deepCR(mask="ACS-WFC-F606W-2-32", inpaint="ACS-WFC-F606W-2-32", device="CPU")
mask, cleaned_image = mdl.clean(image, threshold = thresh)
prob_mask = mdl.clean(image, binary=False)
plt.imshow(cleaned_image)
plt.title('Adjusted Image')
plt.savefig('adjusted.png')
plt.fig()
plt.imshow(image)
plt.title('Original Image')
plt.savefig('orig.png')

