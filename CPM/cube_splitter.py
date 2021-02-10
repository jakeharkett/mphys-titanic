"This script will copy individual images from an image cube"

"FIRST: change the number on line 35 to be the number of the image"
"you want to copy (the number on ds9 minus 1: for example if you"
"want to copy an image that is labeled as image 6 in ds9, type 5"
"in this script), Then run the script:"
"It will ask the user to input the name of the cube, and then"
"will ask for the new name of the output image"

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import warnings

warnings.simplefilter("ignore")

filename = input("Hello!!!\nPlease input the name of your cube (in the form .fits):")

# Sanity check to ensure the image inputted actually exists.
# If not, it will continue to loop until the correct name is inputted
while (not os.path.isfile(filename)):
    print("Cube '%s' not found" % filename)
    filename = input('Please enter the name of the cube (of the form .fits):')

new_filename = input("Please input the name you wish for the new image (in the form .fits):")

hdulist = pyfits.open(filename)
header = hdulist[0].header
data = hdulist[0].data
data = np.nan_to_num(data)
new_data = data[0]

new_data = data[0]

hdu = pyfits.PrimaryHDU(new_data)
hdu.writeto(new_filename)
