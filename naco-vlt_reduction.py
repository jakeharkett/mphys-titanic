"""
===== VLT-NACO REDUCTION SCRIPT - JAMES DAVIES =====
"""

"""
STEPS FOR CODE TO TAKE:
1. Read in images
2. Sequential subtraction
3. Flatfield
4. Crop
5. Combine into cube
6. Standard projectscript data reduction
"""

"""
LOG:
11/11/19 - Start
"""

"""
IMPORT NECESSARY PACKAGES:
"""
import matplotlib.pyplot as plt
import numpy as np
import vip_hci as vip
from vip_hci.preproc import cube_recenter_2dfit, cube_recenter_dft_upsampling, cosmetics, badpixremoval, recentering
from vip_hci.negfc import firstguess
from vip_hci.negfc import cube_planet_free
from vip_hci.negfc import show_corner_plot, show_walk_plot
from vip_hci.negfc import mcmc_negfc_sampling

# VIP 0.9.11 does not use pp_subplots anymore. Use matplotlib instead.
#plots = vip.var.pp_subplots
import hciplot

import math
import pandas as pd
from math import cos
from math import sin
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import figure
import warnings
warnings.simplefilter("ignore")

from vip_hci.metrics.contrcurve import *

from vip_hci.preproc import cosmetics
"""
End of package imports 
"""

"""
============================	INPUT VARIABLES		=================================
"""
#psf_xy = [95,111]			#Centre of psf.
psf_xy = [448,468]			#Centre of psf. #JAMES
reference_xy = [0.1, 0.1]	#Reference bright point defined as a decimal so it loops through the decimal check.
star_xy = [0.1, 0.1] 		#Centre of the star co-ordinates, defined as a decimal so it loops through the decimal check.
averageflux = 0				#Initialises the background flux.				#Just a range of flux's for the first guess algorithm .
flvlup = 15000
flvldown = 1000				#to loop over.

#Initialises the loop counter, loop counter is used to prevent error warning on first trial
#when looping through integer check functions.
Loop_count = 0				
#pxscale_keck = 0.00953 		#Pixel scale for keck 2.
pxscale_keck = 0.02719 		#Pixel scale for NACO L27 camera stated in VIP tutorial
#seperation = 60 			#Radial seperation between synthetic planets in terms of pixels.
seperation = 30 			#Radial seperation between synthetic planets in terms of pixels - JAMES - reduced separation due to the new cropped size
starphot = 0				#Initialises the Starphot parametre.
sigma = 5					#Sets the value of sigma for the contrast curves.
sub_count = 0
im = 0.1
""" End of variables """ 

"""
============================	FUNCTION DEFINITIONS		=================================
"""


"""
Reading fits image names from text file
"""
def ReadImageFilenames(name_input):
	
	#Open file containing filenames and saving the contents
	f = open('{name}_filenames.txt'.format(name = name_input))
	file_names = f.readlines()
	f.close()
	
	# Reads filenames as "xxxx.fits\n" so delete 'newline'
	for i in range(len(file_names)):
		file_names[i] = file_names[i].rstrip()

	return file_names




"""
Starphot calculation calculates the flux from the psf, reads the number of
coadds and the integration time from the headers. The function then uses these
values to calculate the parametre starphot the flux from the images in the cube.
"""

def StarphotCalculation(Number_images, psf_xy, starphot):
			
	for i in range(0, Number_images):	
		hdulist = fits.open(file_names[i], ignore_missing_end = True, verbose = False)
			#Creates list of open files
		Coadd_cube = hdulist[0].header['ESO TPL NEXP'] #JAMES - CHANGE ME #was COADDS
			#Reads the number of coadds from the first fits file
		int_time_cube = hdulist[0].header['EXPTIME'] #ITIME = EXPTIME
			#Reads the integration time for each coadd from the first fits file
	
	hdulist2 = fits.open(initialpsf, ignore_missing_end = True, verbose = False)
	
	int_time_psf = hdulist2[0].header['EXPTIME']
		#Reads the integration time for each coadd from the psf fits file
	Coadd_psf = hdulist2[0].header['ESO TPL NEXP']
		#Reads the number of coadds from the psf fits file
	
	ycord = np.full((Number_images), psf_xy[1], dtype=int)
	xcord = np.full((Number_images), psf_xy[0], dtype=int)
	
	flux = vip.metrics.contrcurve.aperture_flux(
			psf, ycord, xcord, fwhm, ap_factor=1,
			mean=False, verbose=False
			)
		#Returns the sum of pixel values in a circular aperture centered on the input coordinates
	
	psfflux = flux[1]
		#Get second member of 'flux'. Why??
		
	starphot = psfflux * ((int_time_cube * Coadd_cube)/(int_time_psf * Coadd_psf))
	
	return starphot


"""
Readangles reads in the filenames and their parralactic 
angles from their headers and writes them to a text file.
"""

def Readangles(Number_images,name_input):
	
	file_names = ReadImageFilenames(name_input)
	#Loads the filenames into python from the textfile containing them.
	
	Number_images = len(file_names)			#Counts the length of the list.
	angles = np.zeros((1, Number_images))	#Initialises an array for the angles.
	
	#Loop creates a list of open files and reads the 'PARANG' parametre from the headers.
	for i in range(0, Number_images):
		hdulist = fits.open(file_names[i], ignore_missing_end = True, verbose = False)
		#angles[0,i] = hdulist[0].header['ESO TEL PARANG START'] #JAMES EDIT- PARANG CHANGE
		angles[0,i] = (hdulist[0].header['ESO TEL PARANG START'] + hdulist[0].header['ESO TEL PARANG END'])/2 #JAMES EDIT- PARANG CHANGE
	
	#Alters the angles from a 2D array into a 1D list for ease of use.
	
	np.savetxt('{name}_angles.txt'.format(name = name_input), angles)
	# Saves angles into txt
	angles = np.loadtxt('{name}_angles.txt'.format(name=name_input))
	#Loads the angles back into python.

	#JAMES - ADD A NEW FILE FOR THE (N-2)CUBE
	anglez = np.zeros((1, Number_images-2))	#Initialises an array for the angles.
	#Loop creates a list of open files and reads the 'PARANG' parametre from the headers.
	for i in range(0, Number_images - 2):
		hdulist = fits.open(file_names[i], ignore_missing_end = True, verbose = False)
		#anglez[0,i] = hdulist[0].header['ESO TEL PARANG START'] #JAMES EDIT- PARANG CHANGE
		anglez[0,i] = (hdulist[0].header['ESO TEL PARANG START'] + hdulist[0].header['ESO TEL PARANG END'])/2 #JAMES EDIT- PARANG CHANGE
	
	#Alters the angles from a 2D array into a 1D list for ease of use.
	
	np.savetxt('{name}_anglez.txt'.format(name = name_input), anglez)
	# Saves angles into txt
	anglez = np.loadtxt('{name}_anglez.txt'.format(name=name_input))
	#Loads the angles back into python.



"""	
Buildcube initialises an empty array and loads the images into 
this array, removing the badpixels and also subtracting the flatfield.
"""	
#Now number of images -2 due to sequential subtraction
def Buildcube(Number_images, name_input):
	
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names) - 2
	#Counts the length of the list.
	
	#Initialises two cubes to loop through.
	cube0 = np.zeros((Number_images, 1024, 1024))				
	cube = np.zeros((Number_images, 1024, 1024))

	#Initialise the (n+2) cube?
	cube1 = np.zeros((Number_images, 1024, 1024))
	jamescube = np.zeros((Number_images, 100, 100))
	
	flatfield0 = './flat.fits'	#Loads in the flat-field.	
		
	flatfield = vip.fits.open_fits(
		flatfield0, n=0, header=False, 
		ignore_missing_end=True, verbose=True)
	#Opens the flat_field.
	
	Images_loaded = 0	#Initialises a counter.
	
	#Loop opens the images and loads them into the data cube format required.
	#Loop also divides by the flatfield.
	y = 0 #Jitter loop counter
	for i in range (0, Number_images):
	
		Images_loaded = Images_loaded + 1
			#Counter.
		print( "Processing {:d}/{:d}.".format(Images_loaded, Number_images) )
		
		cube0 = vip.fits.open_fits(
			file_names[i], n=0, header=False,
			ignore_missing_end=True, verbose=False)
			#Opens the fits file using VIP.

		cube1 = vip.fits.open_fits(
			file_names[i+2], n=0, header=False,
			ignore_missing_end=True, verbose=False)
			#Opens the (n+2) fits file using VIP.

		#JAMES - DETERMINE CENTER
		if i%10==0 and i!=0: #Everytime the jitter loop resets add 1 to y
			y = y + 1
		z = math.trunc((i-(10*y))/2) + 1 #Gives position number

		if z == 1:
			center = (292,727) #WAS 295,726

		if z == 2:
			center = (292,288) #WAS 295,288

		if z == 3:
			center = (732,288) #WAS 732,288

		if z == 4:
			center = (732,726) #WAS 732,726

		if z == 5:
			center = (513,507) #WAS 513,507
		
		for j in range(0, 1024):
			for k in range(0, 1024):
				cube[i][j][k] = (cube0[j][k] - cube1[j][k]) / flatfield[j][k] #sequential subtraction and flatfielding
	
		jamescube[i] = cosmetics.frame_crop(array = cube[i], size = 100, cenxy=center, force=False, verbose=True)
		
		if Images_loaded == Number_images:
			print( "\nImages loaded in.\n" )

	#Removes the 'bad pixels' in the image using a VIP function to do so.
	print( "Removing the badpixels in the Cube...\n" )
	
	cube = vip.preproc.badpixremoval.cube_fix_badpix_isolated(
		jamescube, bpm_mask=None, sigma_clip=3,
		num_neig=5, size=5, protect_mask=False,
		radius=30, verbose=True, debug=False)		

	#Writes the created cube to a fits file.	
	vip.fits.write_fits('Newcube.fits', jamescube,  verbose=True)	
	
	#Loads the angles and appends them to the cube, then overwrites the previous cube
	#with the appended one.
	   
	angles = np.loadtxt('{name}_angles.txt'.format(name=name_input))
	vip.fits.info_fits('Newcube.fits')
	
		


"""
Removeimage removes the bad images by erasing 
the name of the file from the list of filenames 
"""

def Removeimage(Number_images,name_input):

	#Loads in the filenames.
	file_names = ReadImageFilenames(name_input)
	user_input = input(
		"Type the file numbers (eg. 1,4,5.. ) of the images you would like to remove\n Seperated by single commas\n")
	
	Bad_images = [int(i) for i in user_input.split(',') if i.isdigit()]
	Bad_images = sorted(Bad_images, reverse =True)	#Sorts the numbers in descending order.
	
	#Loop removes the filenames corresponding to the numbers entered from the list of filenames.
	for i in range(0,len(Bad_images)):
		file_names = np.delete(file_names, (Bad_images[i]-1))
	
	return file_names

"""		
Checkfunction allows for yes no answers and prompts the 
user to re-enter if initially enetered incorrectly.

Yes = 0
No = 1

"""

def Checkfunction():

	check=2		
	while check == 2:	#Ensures that this will be looped over unless a yes or no is entered.
	
		yes=set(['yes', 'y', 'Y','Yes','YES'])	#Sets any of these entries to equal yes.
		no=set(['no', 'n', 'N','No','NO'])	#Sets any of these entries to equal no.
		exit=set(['EXIT','exit','exit()','Exit()','Exit','EXIT()'])
		question=input('Yes (y) or no (n)?\n')
		if question in yes:
			check=0
			
		if question in no:
			check=1	
		
		if question in exit:
			raise SystemExit()
		
		#An error if loop to prompt the user to change their response.	
		if question not in yes and question not in no:
			print( '\nUnrecognised answer\n' )
		
	return check

"""
Contrastcurvedata writes the contrast curve results into 
a text file, saves and then loads it back, plots the LLSG
and the PCA data on the same graph and saves it to a file
"""
"""
Note: This function is terribly inefficient, could be much 
better looped over but haven't had the time to do so.
"""
def Contrastcurvedata(PCA_contrast_curve, LLSG_contrast_curve, pxscale_keck, sigma,name_input):
	
	#Saves the PCA,ADI and LLSG curve outputs and reads them back in. 
	np.savetxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name_input),PCA_contrast_curve)
	np.savetxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name_input),LLSG_contrast_curve)
	

	#Loads the PCA,ADI and LLSG curve outputs and reads them back in. 
	PCA_data = np.loadtxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name_input))
	LLSG_data = np.loadtxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name_input))


	#JAMES - CHANGE SIZE FROM 500 TO 150 AND MAKE IT VARIABLE
	jamessize = 38 #max size i can make this
	#jamessize = 100 #should be 150 i think?
	#jamessize = 20
	
	#Initialises the PCA variables (500 being the size).
	PCA_Distance = np.zeros(jamessize)
	PCA_Sensitivity_Gauss = np.zeros(jamessize)
	PCA_Sensitivity_Student = np.zeros(jamessize)
	LLSG_Distance = np.zeros(jamessize)
	LLSG_Sensitivity_Gauss = np.zeros(jamessize)
	LLSG_Sensitivity_Student = np.zeros(jamessize)

	"""
	PANDAS DATAFRAME: (IF STUDENT == TRUE)
	0 - sensitivity_gaussian
	1 - sensitivity_student
	2 - throughput
	3 - distance
	4 - distance_arcsec
	5 - noise
	6 - sigma corr

	PANDAS DATAFRAME: (IF STUDENT == FALSE)
	0 - sensitivity_gaussian
	1 - throughput
	2 - distance
	3 - distance_arcsec
	4 - noise
	"""

	#JAMES - NEW FUNCTION - Updated pandas database columns
	#Loop loads the PCA outputs into the PCA variables.
	for a in range (0,jamessize):
		PCA_Distance[a] = PCA_data[a,3] * pxscale_keck
		PCA_Sensitivity_Gauss[a] = PCA_data[a,0]
		PCA_Sensitivity_Student[a] = PCA_data[a,1]	

		#Loop loads the LLSG outputs into the LLSG variables.
		LLSG_Distance[a] = LLSG_data[a,3] * pxscale_keck
		LLSG_Sensitivity_Gauss[a] = LLSG_data[a,0]
		LLSG_Sensitivity_Student[a] = LLSG_data[a,1]

	#OLD FUNCTION - values for pandas database have changed since this was written
	"""
	#Loop loads the PCA outputs into the PCA variables.
	for a in range (0,jamessize):
		PCA_Distance[a] = PCA_data[a,0] * pxscale_keck
		PCA_Sensitivity_Gauss[a] = PCA_data[a,2]
		PCA_Sensitivity_Student[a] = PCA_data[a,3]	

		#Loop loads the LLSG outputs into the LLSG variables.
		LLSG_Distance[a] = LLSG_data[a,0] * pxscale_keck
		LLSG_Sensitivity_Gauss[a] = LLSG_data[a,2]
		LLSG_Sensitivity_Student[a] = LLSG_data[a,3]
	"""
	

	#Plotting all 3 on one curve: 
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#JAMES - deleting mysterious cross
	"""
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	"""
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[PCA_Legend,LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.

	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	#plt.ylim((0.000001,0.01))#JAMES - LIMIT NEEDS TO BE CHANGED FOR L-BAND
	
	#Creates the variables that will be plotted.
	Curve = [0,0]	#Initialises the Curve variable.
	Curve[1] = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	Curve[0] = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	plt.savefig('new_Contrast_curves_no_ADI_{star_name}.png'.format(star_name=name_input))

	plt.show(Curve)
	
	#LLSG contrast curve 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#JAMES - deleting mysterious cross
	"""
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	"""
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	#plt.ylim((0.000001,0.01))#JAMES - LIMIT NEEDS TO BE CHANGED FOR L-BAND
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	#Saves the figure and then shows the plot. 
	plt.savefig('new_LLSG_curves_{star_name}.png'.format(star_name=name_input))

	plt.show(Curve)
	
	#PCA contrast curve. 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#JAMES - deleting mysterious cross
	"""
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	"""
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	plt.legend(handles=[PCA_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	#plt.ylim((0.000001,0.01)) #JAMES - LIMIT NEEDS TO BE CHANGED FOR L-BAND
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	
	#Saves the figure and then shows the plot. 
	plt.savefig('new_PCA_curve_{star_name}.png'.format(star_name=name_input))
	plt.show(Curve)


"""
Function takes the cube and the PSF (plus its centre)
and fits a 2d Gaussian into it.
Using that, it recenters the cube and plots the result.
Returns the re-centered cube, and the x/y shifts
"""

def Gaussian_2d_Fit( psf, cube_orig, centre):

	#OLD FUNCTION
	"""
	gauss = vip.var.fit_2dgaussian(
		psf, crop=True, cropsize=30, cent=(psf_xy[0], psf_xy[1]),
		full_output=True, debug=True
		)
	"""

	#JAMES -NEW FUNCTION
	gauss = vip.var.fit_2d.fit_2dgaussian(
		psf, crop=True, cropsize=30, cent=(psf_xy[0], psf_xy[1]),
		full_output=True, debug=True
		)
	
	#Calculates the FWHM from the VIP gaussian function.
	
	print( gauss )
	
	# 'gauss' is a pandas dataframe
	# Columns are: 
	#	0-centroid_y		1-centroid_x	2-fwhm_y	3-fwhm_x	4-amplitude		5-theta
	# Hence take columns 2 and 3 for FWHM y and x components
	fwhm_x = gauss.iloc[0,3]
	fwhm_y = gauss.iloc[0,2]
	
	
	fwhm = np.mean([fwhm_x, fwhm_y])
	print( "FWHM = ", fwhm )

	fwhm_int = math.ceil(fwhm)
	fwhm_int = int(fwhm_int)
	print( "Round up = ", fwhm_int )
	
	#Shy1 = Shift in y, shx1 = shift in x. Calculates the shifts needed to centre the image.
	#Negative = True (Bright spot at the centre) Negative = False (No bright spot at centre).
	cube1, shy1, shx1 = cube_recenter_2dfit(
		cube_orig, (centre[0], centre[1]), fwhm,
		model='gauss', nproc = 1, subi_size = (fwhm_int), #JAMES EDIT TO 3*
		negative=False, full_output=True, debug=False
		)
		
	#Plots the shifts in x and y used to centre the images.
	plt.figure(figsize = (8,3))								
	plt.plot(shy1, 'o-', label='shifts in y', alpha=0.5)	
	plt.plot(shx1, 'o-', label='shifts in x', alpha=0.5)
	plt.legend(loc = 'best')
	plt.grid('on')
	plt.suptitle('Positive 2d Gaussian fit')
	plt.show()
	
	return cube1, shy1, shx1, fwhm
	
		
""" 
End of defining functions 
"""


""" 
************************************************************************************************
=============================	CODE STARTS HERE =============================
************************************************************************************************
"""

print("\n\n\n ================= Project Script - VLT-NACO - JAMES DAVIES =================")
print("\n Version 11_py3 - VLT-NACO - JAMES DAVIES ADAPTATION\n")
print("VIP version is:",vip.__version__)
print("\n\n")
	
"""
------ 1 - GET STAR NAME
"""
	
name_input = input('Name of the star (example = HIP_544):  ')
Centring_loop = 0 # Initialise loop to allow for re centring

while Centring_loop == 0:
	
	#Loads the filenames into python and counts the number of files.
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names)
	
	
	"""
	------ 2 - CREATE DATA CUBE
	"""
	build_check = 0
	#Initialises the variable which is used to build the cube (or not).
	print( "Would you like to create the data cube?")
	build_check = Checkfunction()

	#Loop runs the Readangles and Buildcube functions, building the data cube.
	if build_check == 0:
		file_names = ReadImageFilenames(name_input)
		Number_images = len(file_names)
		Readangles(Number_images,name_input)
		Buildcube(Number_images,name_input)

	
	print( "Open DS9 and look through the created data cube at each frame, taking note of the bad frames" )
	remove_check = 0
	
	
	"""
	------ 3 - REMOVE BAD IMAGES
	"""
	print( "Would you like to remove any image?" )
	remove_check = Checkfunction()

	#Loop runs the Removeimage function and then overwrites the new shortened list to
	#the original text file. 
	if remove_check == 0:
		
		file_names = Removeimage(Number_images,name_input)
		file_names = file_names.tofile('{name}_filenames.txt'.format(name=name_input), sep ="\n", format='%s')
		file_names = ReadImageFilenames(name_input)
		Number_images = len(file_names)
		Readangles(Number_images,name_input)
		Buildcube(Number_images,name_input)
	

	"""
	------ 4 - RECENTER LOOP #JAMES - FIX ME FOR NEW CROPPED SIZE
	"""
	#Loads the psf and the cube fits into python.
	initialpsf = './psf.fits'			
	cube = './Newcube.fits' 			


	Loop_counttwo = 1
	# Need to re-set star_xy so it asks for coordinates again 
	# and need to stop it asking to center loop again and go straight to getting coordinates and exiting!
	
	Loop_countthree = 0
	#auto_find_coords = 0
	# Needed to make sure loop does not allow you to keep running loop over centre coordinates

	while Loop_counttwo > 0:
	
		star_xy = [0.1, 0.1]
		
		print("Using DS9, open any image in Newcube.fits")
		star_xy[0] = input("Input the x co-ordinate of the central position of the star (integer number): ")
		star_xy[0] = int(star_xy[0])
		
		star_xy[1] = input("Input the y co-ordinate of the central position of the star (integer number): ") #JAMES-EDIT FIXED TYPO
		star_xy[1] = int(star_xy[1]) 
		
		
		#Opens the cube (using VIP) with cube_orig as HDU:0 and calls it cube_orig, 
		#and the parallactic angles from a text file.
		#Uses VIP to also open the point spread function previously loaded in.
		
		cube_orig = vip.fits.open_fits(cube)
		angs = np.loadtxt( '{name}_angles.txt'.format(name = name_input) )
		psf = vip.fits.open_fits(initialpsf, n=0, header=False, ignore_missing_end=True, verbose=True)
	
		print( "Fitting a 2D gaussian to centre the images..." )
		cube1, shy1, shx1, fwhm = Gaussian_2d_Fit( psf, cube_orig, star_xy)
		print("FWHM = ",fwhm) #JAMES-CHECK
		
		#Uses VIP's 2D gaussian fitting algorithm to centre the cube.
		
		

		# Created a loop to get over the annoying problem of choosing a perfect centre
		# Will look through various centres and decides if it's a good fit or not.
		# Maybe use its own function??
		
		"""
		Loop through new coordinates ========================================================================
		"""
		
		
		
		if Loop_countthree == 0:  #only allows the loop over centre coordinates on first iteration of control loop
			
			print( "If you are having difficulty finding the best centre point," )
			print( "then you can loop through a few co-ordinates close to the original value" )
			print( "to find a better centre" )
			check_gauss = 0
			print( "Would you like to loop over a region close to the centre value?" )
			check_gauss = Checkfunction()
			
			# If you want to loop with sightly different coordinates:
			if check_gauss == 0:
				Loop_centre = [0,0]
				Loop_centre[0] = star_xy[0] - 2
				Loop_centre[1] = star_xy[1] - 2
				Loop_break=1 # allows for the loop to be broken to save time if correct values found quickly
				print( "New center x =",Loop_centre[0] )
				print( "New center y =",Loop_centre[1] )
				for Loop_centre[0] in range (star_xy[0] -2, star_xy[0] +2):
	
					for Loop_centre[1] in range(star_xy[1]-2, star_xy[1] +2):
		
						print( "Using the value [{:f},{:f}] as the centre".format(Loop_centre[0], Loop_centre[1]) )
						
						cube1, shy1, shx1, fwhm = Gaussian_2d_Fit( psf, cube_orig, Loop_centre)
		
						sumshy1 = 0
						sumshx1 = 0
						for i in range (0, Number_images):
							sumshy1 = sumshy1 + shy1[i]
							sumshx1 = sumshx1 + shx1[i]
			
						averageshy = sumshy1/Number_images
						averageshx = sumshx1/Number_images
						shifty = np.zeros((Number_images))
						shiftx = np.zeros((Number_images))
						Bad_count = 0
		
						Loop_centre[1] = Loop_centre[1] + 1
		
						for i in range (0, Number_images):
							shifty[i] = shy1[i] - averageshy
							shiftx[i] = shx1[i] - averageshx
		
							if shiftx[i] > 5 or shifty [i] > 5:
								Bad_count = Bad_count + 1
			
						if Bad_count > Number_images/5:
							print( "This centre is not very good for VIP's Gauss fit" )
		
						if Bad_count < Number_images/5:
							print( "This centre is a good fit" )
				
						print( "Exit loop?" )
						Loop_break = Checkfunction()
						if Loop_break == 0:
							break
				
					Loop_centre[0] = Loop_centre[0] + 1
					if Loop_break == 0:
						break
						
			Loop_counttwo = Loop_counttwo + 1
		"""
		End of coord loop ========================================================================
		"""
		
		"""
		# If 
		if check_gauss == 1:
			cube1, shy1, shx1, fwhm = Gaussian_2d_Fit( psf, cube_orig, star_xy)
			Loop_counttwo = 0
		"""
		
		Loop_counttwo = 0
		Loop_countthree = Loop_countthree + 1
		
	
	#Writes the values of the centered cube into a fits file.
	vip.fits.write_fits('centeredcube_{name}.fits'.format(name=name_input), cube1, verbose=True)	

	cube = cube1
	#Loads up the centered cube.
	#Plots the original cube vs the new centered cube.

	#JAMES - CROP
	im1 = vip.preproc.cosmetics.frame_crop(cube_orig[0], 30, verbose=True) #JAMES EDIT TO 30
	im2 = vip.preproc.cosmetics.frame_crop(cube[0], 30, verbose=True)
	
	#Plotting the original frame vs the new recentered frame
	hciplot.plot_frames( 
		(im1, im2), 
		label = ('Original first frame', 'First frame after recentering'), 
		grid = True, 
		size_factor = 4
		)
	
	print( "Open DS9 and look through the centred data cube at each frame making sure it is centred.")
	print( "If you're not happy with it, redo centring" )
	print( "Redo Centring?" )
	Centring_loop = Checkfunction()


"""
------ 5 - REFERENCE COORDINATES
"""
	
print( "Use DS9 to obtain the reference co-ordinates\n" )

reference_xy[0] = input("Input the x co-ordinate of the reference/candidate position (Integer number): ")
reference_xy[0] = int(reference_xy[0])

reference_xy[1] = input("Input the y co-ordinate of the reference/candidate position (Integer number): ")
reference_xy[1] = int(reference_xy[1])

"""
----5.5 - JAMES RESET STUPID GAUSSIAN >:(
"""
cube = cube_orig


"""
------ 6 - PCA
"""
print( 'Will now reduce the data using full frame PCA' )

#Reduce the image using full frame PCA then plots the image..
# Optimal number of PCs are calculated

anglez = np.loadtxt( '{name}_anglez.txt'.format(name = name_input) ) #new list of (n-2) angles

ffpca_output = vip.pca.pca_fullfr.pca(
	cube, anglez, fwhm = fwhm,
	source_xy = (reference_xy[0], reference_xy[1]),
	mask_center_px = None,  ncomp = (1, len(cube) ,1), 
	full_output=True
	)
	
# vip_hci.pca.pca_fullfr.pca returns this outputs:
# (if full_output=True and source_xy != None)
# 0 - Residuals cube (for STIM map input)
# 1 - PCA Frame
# 2 - Pandas dataframe with PCs, S/Ns and fluxes
#
#<class 'numpy.ndarray'> 10,1024,1024
#<class 'numpy.ndarray'> 1024,1024
#<class 'pandas.core.frame.DataFrame'> 10


fr_pca1 = ffpca_output[1]
resid_cube_derot = ffpca_output[0]

		

hciplot.plot_frames( 
		fr_pca1, 
		label = 'FFPCA reduction of {name}'.format(name=name_input), 
		grid = False, 
		size_factor = 5
		)

#Loop asks user if they would like to save the image.
print( 'Would you like to save the PCA reduced images?' )
questionsave = Checkfunction()

if questionsave == 0:	
	vip.fits.write_fits('new_FFPCA_{name}.fits'.format(name=name_input), fr_pca1, verbose=True) 


"""
------ 7 - LLSG
"""
llsgloop = 0 #Initialises the LLSG loop variable.

#Loop to allow ranks of llsg to be altered if the image is not optimal.
while llsgloop ==0:	
	
	Rank = input("What rank of llsg would you like to do? (default = 6)\n")
	Rank = int(Rank)
	
	#Reduces the image using the LLSG algorithm then plots the image.
	print(fwhm)
	fr_llsg = vip.llsg.llsg(cube, anglez, fwhm=fwhm, rank = Rank, thresh = 2, max_iter = 20, random_seed = 10)
	hciplot.plot_frames(fr_llsg, label ='LLSG reduced image of {name}'.format(name=name_input))
	
	
	#Loop asks user if the would like to save the image.
	print( 'Would you like to save the LLSG reduced image?\n' )
	questionsave = Checkfunction()
	
	if questionsave == 0:
		vip.fits.write_fits('new_LLSG_{name}.fits'.format(name=name_input), fr_llsg, verbose=True) 
	
	#Checks if the user would like to use a different rank and stays in while loop if so.	
	print( "Would you like to try a different value of rank?\n" )
	llsgloop = Checkfunction()


"""
------ 8 - STARPHOT CALCULATION
"""
#Starphot is the  total aperture flux from the star during
#the integration time (if no coronograph).
starphot = StarphotCalculation(Number_images-2, psf_xy, starphot)
#starphot_array = np.zeros(500)
starphot_array = np.zeros(150) #JAMES JAMES-EDIT

#for i in range (0,500):
for i in range (0,150):
	starphot_array [i] = starphot

np.savetxt('{name}_starphot'.format(name=name_input), starphot_array)


"""
------ 10 - CONTRAST CURVES
"""
#JAMES - FIX ME!!!!
#psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=101, verbose=True, debug=True, threshold=None, mask_core=None,)
psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=31, verbose=True, debug=True, threshold=None, mask_core=None,)
#psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=30, threshold=None, mask_core=None)

fwhm = fwhm * 2 #JAMES - FOR SOME REASON THIS NEEDS TO BE LARGER??	

#JAMES -EDIT TO MAKE STUDENT TRUE?
PCA_contrast_curve = vip.metrics.contrcurve.contrast_curve(
	cube, anglez, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=1, algo=vip.pca.pca, 
	debug=True,save_plot='PCA', verbose=2, student=True,
	ncomp=40 #optimal number is 10? not 40?
	)

LLSG_contrast_curve = vip.metrics.contrcurve.contrast_curve(
	cube, anglez, psf, fwhm,pxscale_keck, starphot,
	sigma=sigma, nbranch=1, algo=vip.llsg.llsg, student=True,
	debug=True,save_plot='LLSG',verbose=2
	)

Contrastcurvedata(PCA_contrast_curve, LLSG_contrast_curve, pxscale_keck, sigma, name_input)



"""
------ 11 - STIM MAPS #JAMES - JORGE BROKE ME :'(
"""
"""
print("Computing STIM map...")
stim_map = vip.metrics.compute_stim_map(resid_cube_derot)


hciplot.plot_frames( 
		stim_map, 
		label = 'STIM map of {name}'.format(name=name_input), 
		grid = False, 
		size_factor = 5
		)

#Loop asks user if they would like to save the image.
print( 'Would you like to save the image?' )
questionsave = Checkfunction()

if questionsave == 0:	
	vip.fits.write_fits('STIM_map_{name}.fits'.format(name=name_input), stim_map, verbose=True)
	
"""
"""
------ 12 - ANDROMEDA MAP
"""






print( "===========================	End of the program	===========================" )





# Some space here