'''
This program measures the background of the image and looks for
objects that are 50X brighter than this background with a
significance of 5sig. Be warned that this program has a tendancy
to pick up bright speckles as well as planets, although the planet
should be by far the most prominent signal. Use ds9 to find the
approx position of the planet and then this program to find a more
accurate position. Ignore the other data such as flux and magnitude
that appears in the table as this is not calibrated.
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import math
import pandas as pd
import sys
import argparse

from astropy.io import fits
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.visualization import astropy_mpl_style

from photutils.background import MADStdBackgroundRMS
from photutils.psf import IntegratedGaussianPRF
from photutils.psf import DAOPhotPSFPhotometry
from photutils import DAOStarFinder
from photutils import CircularAperture
from photutils import aperture_photometry
from termcolor import colored
import warnings

warnings.simplefilter("ignore")

print('Good Morning!\n')
image_file = input('Please enter the name of the image (of the form .fits):')

# Sanity check to ensure the image inputted actually exists.
# If not, it will continue to loop until the correct name is inputted
while (not os.path.isfile(image_file)):
    print(colored("Image '%s' not found" % image_file,"red"))
    image_file = input('Please enter the name of the image (of the form .fits):')


# Input parameters for the detection code
signal_noise = input('Please enter the signal to noise threshold you would like:')
signal_noise = float(signal_noise)

sigma_psf = input("Please enter the value of sigma you wish the psf to have:")
sigma_psf = float(sigma_psf)

while(signal_noise <= 0 ):
    print(colored('Signal to noise threshold must be greater than 0','red'))
    signal_noise = input('Please enter the signal to noise threshold you would like:')
    signal_noise = float(signal_noise)

while(sigma_psf <= 0):
    print(colored('Sigma must be greater than 0','red'))
    sigma_psf = input("Please enter the value of sigma you wish the psf to have:")
    sigma_psf = float(sigma_psf)

hdulist = fits.open(image_file, ignore_missing_end=True)
hdu_data = hdulist[0].data
hdulist.close()

bkgrms = MADStdBackgroundRMS()
std = bkgrms(hdu_data)
thresh = signal_noise * std
fwhm_sigma = sigma_psf * gaussian_sigma_to_fwhm
fitshape = int(3 * np.ceil(fwhm_sigma) // 2 * 2 + 1)

stfind = DAOStarFinder(threshold=thresh, fwhm=fwhm_sigma)
sources = stfind(hdu_data)
for col in sources.colnames:
    sources[col].info.format = '%.8g'
print(sources)

plt.style.use(astropy_mpl_style)
image_data = fits.getdata(image_file, ext=0, ignore_missing_end=True)

plt.figure()
plt.imshow(image_data, cmap='gray')
plt.gca().invert_yaxis()

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=16.)
apertures.plot(color='red', lw=1.5)


print(colored('Ensure you have used ds9 to check the rough positions and identify the planet\n','red'))

print('==================== End of the program =======================')
