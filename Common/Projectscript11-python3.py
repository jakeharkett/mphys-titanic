""" 
**************************************NOTES****************************************

Script last updated 29/05/19
this script works for VIP 0.9.8 in Python2 environment (current version on analysis)

This is different from Projectscript9.py

**************************************INSTRUCTIONS**********************************

Enter the name of the star to run the script, e.g. HIP_99542. 
Put this folder in the epoch file. 

Please upload contrcurve.py, cosmetics.py, fakecomp.py and shapes.py when running this script (avaliable in DangerZone)

**************************************UPDATES**************************************

01/06/17
	1. Option to break from loop when looping over different centre postions without
	   having to wait till all positions tried
	2. Stops users enter more six planets at planet injection stage, recomends use of 6
	
05/06/17
	1. Includes allowing for coordinate re-entry after looping over centre positions
	2. allow user to enter what multiple of background flux is used for planet brightness

06/06/17
	1. Clearer instructions and error checking when entering planet flux
	2. Background flux check now takes four 100x100 pixel grids from each image in the 
	   data cube. If the average flux in any of each of this grids is greater than 1 stdev
	   away from the population mean that grid is rejected and not used to build a mean
	   background flux
	3. Minor adjustment to Checkfuntion() allowing for more input options
	4. Changed Checkfunction to allow for exit to be called.	

07/06/17
	1. Addition of loop to allow planet injection phase to run mutliple times.
	2. Loop to redo entirly the re-centring if ds9 inspection shows problem
	3. Made if more obvious when you need to enter reference coordinates
	4. fixed bug in coordinate re-entry
	5. Made Image removal function more efficient and removed the necessity to know how many
	   images were going to be removed before entering numbers.
	6. Background flux now takes grids within 0.5 stdev of mean
08/06/17
	1. Removed Planet injection loop due to bug
	2. The .png files of each reduction were outputing blank files. No way was found to 
	   fix this therefore the code was removed. still saved as .fits but no .png.
	   This appears to be beacuse it is saving wrong figure but this is due to a VIP
	   hence no solution was found.
14/06/17
	1. Fixed the ininate loop when asking to re enter coodinates
	2. Fixed the re entered coordinates not being correctly over writtten

28/05/18
	1. Update VIP from 0.5.4 to 0.9.4
	2. Note that any script older than Projectscript 7 will not work now as structure of packages are different 
	3. Update the script due to the difference in packages

28/06/18 
	1. Add planet subtraction as a loop
	2. Edited the script hence user can run the script by entering the name

28/05/19 
	1. Avoid the crash of contrast curve by doing : 
		import fakecomp.py cosmetics.py shapes.py contrcurve.py (the script I have uploaded on analysis)
		from contrcurve import *
		change "vip.metrics.contrcurve.contrast_curve" to "contrast_curve"

29/05/19 
	1. Fixed bugs due to the updates
	2. Disable Planet subtraction (there should a bug in matplotlib)


----
02/10/19
	1. Porting to Python3
		Changed 'print' to 'print()'
		Changed 'raw_input()' to 'input()'
	
03/10/19
	1. Updating VIP from 0.9.8 to 0.9.11
		pp_subplots no longer supported. Use HCIplot module
		Previous scripts will no longer work, Python 3 required!!
		Function open_adicube() no longer in VIP. Use fits.open_fits()
		Function vip.var.fit_2dgaussian() returns pandas dataframe, 
			output columns in different order now.
		
	2. Porting to Python 3
		Reading filenames changed from np.loadtxt() to standard Python open()
		because Python 3 treats strings and bytes differently
	
	3. Changed requiring-int-loops to simply reading input and using int()
	
04/10/19
	1. Cleaning code
		backgroundflux function is repeated!
		Adding spaces and indentation
		
	2. Porting to VIP 0.9.11
		Function vip.pca.pca_optimize_snr() doesn't work anymore.
		
07/10/19
	Conflict of names: HIP_number and HIPnumber
	Created new function to fit 2d gaussian
		
		
		
		

**************************************************************************************   
	
"""

"""
Script designed to build a data cube from fits files with the aid of the VIP module
in python, reduce the images using the PCA and LLSG method of imaging reduction. 
Script also then calculates the flux from the exposure times of the fits files and 
then generates contrast curves based on the data cube produced and this value of flux.


STUCTURE

	1. Packages
	
	2. Input variables
	
	3. Functions
""" 


"""
============================	PACKAGES	=================================
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
"""
End of package imports 
"""





"""
============================	INPUT VARIABLES		=================================
"""
psf_xy = [95,111]			#Centre of psf.
reference_xy = [0.1, 0.1]	#Reference bright point defined as a decimal so it loops through the decimal check.
star_xy = [0.1, 0.1] 		#Centre of the star co-ordinates, defined as a decimal so it loops through the decimal check.
averageflux = 0				#Initialises the background flux.				#Just a range of flux's for the first guess algorithm .
flvlup = 15000
flvldown = 1000				#to loop over.

#Initialises the loop counter, loop counter is used to prevent error warning on first trial
#when looping through integer check functions.
Loop_count = 0				
pxscale_keck = 0.00953 		#Pixel scale for keck 2.
seperation = 60 			#Radial seperation between synthetic planets in terms of pixels.
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
		Coadd_cube = hdulist[0].header['COADDS']
			#Reads the number of coadds from the first fits file
		int_time_cube = hdulist[0].header['ITIME']
			#Reads the integration time for each coadd from the first fits file
	
	hdulist2 = fits.open(initialpsf, ignore_missing_end = True, verbose = False)
	
	int_time_psf = hdulist2[0].header['ITIME']
		#Reads the integration time for each coadd from the psf fits file
	Coadd_psf = hdulist2[0].header['COADDS']
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
Backgroundflux uses a 100 by 100 square away from the centre 
of the images to get a rough guess of the background noise. 
It takes four 100x100 squares for each image in the data cube 
at oposing points of the image. It then checks that each cube has an
average flux within 0.5 stdev of the population mean. If this criteria
is failed the square is thown out. The average flux is calcuated from 
the remaining squares
"""

def backgroundflux(cube,averageflux):

	number_images = len(cube)
		#calculates number of image in cube
		
	totalflux = np.zeros(shape = [4, number_images])
		#Initialises the total flux array variable.
		
	averageflux_array = np.zeros(shape = [4, number_images])
		#Initialises an array to contain average of each 100x100 grid
													
	for k in range (0, (number_images)):
		for i in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[0,k] = totalflux[0,k] + cube[k,i,j]
				
		for i in range (100, 201): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[1,k] = totalflux[1,k] + cube[k,i,j]
			
		for i in range (800, 901): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(800, 901):
				totalflux[2,k] = totalflux[2,k] + cube[k,i,j]
			
		for i in range (900, 801): #Loops over a 100 by 100 pixel grid to calculate the total flux in that grid.
			for j in range(100, 201):
				totalflux[3,k] = totalflux[3,k] + cube[k,i,j]
		
		for i in range(0, 4):
			averageflux_array[i,k] = totalflux[i,k]/10000
			#Calculates the average flux in each grid by dividing by number of pixels (10000)
			
	upper_bound = np.mean(averageflux_array) + (np.std(averageflux_array)/2)
	lower_bound = np.mean(averageflux_array) - (np.std(averageflux_array)/2)
		#Creates upper and lower bounds of 0.5 stdev above and below population mean
	
	
			
	check_array = np.ones(shape = [4, number_images])
	#Initialises array of ones so Boolean logic will accept as 'true'. Same size and shape as averageflux_array
	
	for k in range (0, (number_images)):
		for i in range(0, 4):
			if averageflux_array[i,k] > upper_bound or averageflux_array[i,k] < lower_bound: 
				check_array[i,k] = 0
				#If flux out of bounds changes check_array data point to zero so boolean logic would give false
				
	Loop_count = 1 #counts how many images are used
	totalflux = 0 #Initialises total flux variable that counts all flux
	
	for k in range (0, (number_images)): 
		for i in range(0, 4):
			if check_array[i,k]:
				#If check_array point [i,k] is true, it adds the corresponding averageflux_array point and adds to totalflux
				totalflux = totalflux + averageflux_array[i,k]
				Loop_count = Loop_count + 1
			
	averageflux = totalflux/Loop_count
	#Divides the totalflux by the number of images used.
	
	print( "Number of grids used = \n", Loop_count )
	print( "Average flux = \n", averageflux )
	print( "Standard Deviation = \n", np.std(averageflux_array) )
	Loop_count = 0
	
	return averageflux
	
	
	
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
		angles[0,i] = hdulist[0].header['PARANG']
	
	#Alters the angles from a 2D array into a 1D list for ease of use.
	
	np.savetxt('{name}_angles.txt'.format(name = name_input), angles)
	# Saves angles into txt
	angles = np.loadtxt('{name}_angles.txt'.format(name=name_input))
	#Loads the angles back into python.



"""	
Buildcube initialises an empty array and loads the images into 
this array, removing the badpixels and also subtracting the flatfield.
"""	

def Buildcube(Number_images,name_input):
	
	file_names = ReadImageFilenames(name_input)
	Number_images = len(file_names)
	#Counts the length of the list.
	
	#Initialises two cubes to loop through.
	cube0 = np.zeros((Number_images, 1024, 1024))				
	cube = np.zeros((Number_images, 1024, 1024))				
	
	flatfield0 = './flat_Kp.fits'	#Loads in the flat-field.	
		
	flatfield = vip.fits.open_fits(
		flatfield0, n=0, header=False, 
		ignore_missing_end=True, verbose=True)
	#Opens the flat_field.
	
	Images_loaded =0	#Initialises a counter.
	
	#Loop opens the images and loads them into the data cube format required.
	#Loop also divides by the flatfield. 
	for i in range (0, Number_images):
	
		Images_loaded = Images_loaded + 1
			#Counter.
		print( "Processing {:d}/{:d}.".format(Images_loaded, Number_images) )
		
		cube0 = vip.fits.open_fits(
			file_names[i], n=0, header=False,
			ignore_missing_end=True, verbose=False)
			#Opens the fits file using VIP.
		
		for j in range(0, 1024):
			for k in range(0, 1024):
				cube[i][j][k] = cube0[j][k] / flatfield[j][k]
		
		if Images_loaded == Number_images:
			print( "\nImages loaded in.\n" )
	
	#Removes the 'bad pixels' in the image using a VIP function to do so.
	print( "Removing the badpixels in the Cube...\n" )
	
	cube = vip.preproc.badpixremoval.cube_fix_badpix_isolated(
		cube, bpm_mask=None, sigma_clip=3,
		num_neig=5, size=5, protect_mask=False,
		radius=30, verbose=True, debug=False)		
		
	#Writes the created cube to a fits file.	
	vip.fits.write_fits('Newcube.fits', cube,  verbose=True)	
	
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
	

	#Saves the PCA,ADI and LLSG curve outputs and reads them back in. 
	#PCA_data = np.loadtxt('new_PCA_{star_name}_curve_outputs'.format(star_name=name_input))
	#LLSG_data = np.loadtxt('new_LLSG_{star_name}_curve_outputs'.format(star_name=name_input))
	
	#Initialises the PCA variables (500 being the size).
	PCA_Distance = np.zeros(500)
	PCA_Sensitivity_Gauss = np.zeros(500)
	PCA_Sensitivity_Student = np.zeros(500)
	LLSG_Distance = np.zeros(500)
	LLSG_Sensitivity_Gauss = np.zeros(500)
	LLSG_Sensitivity_Student = np.zeros(500)

	#Loop loads the PCA outputs into the PCA variables.
	for a in range (0,500):
		PCA_Distance[a] = PCA_data[a,0] * pxscale_keck
		PCA_Sensitivity_Gauss[a] = PCA_data[a,2]
		PCA_Sensitivity_Student[a] = PCA_data[a,3]	

	#Loop loads the LLSG outputs into the LLSG variables.
		LLSG_Distance[a] = LLSG_data[a,0] * pxscale_keck
		LLSG_Sensitivity_Gauss[a] = LLSG_data[a,2]
		LLSG_Sensitivity_Student[a] = LLSG_data[a,3]
	
	#Plotting all 3 on one curve: 
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[PCA_Legend,LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.

	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = [0,0]	#Initialises the Curve variable.
	Curve[1] = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	Curve[0] = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	savefig('new_Contrast_curves_no_ADI_{star_name}.png'.format(star_name=name_input))

	plt.show(Curve)
	
	#LLSG contrast curve 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	LLSG_Legend = mpatches.Patch(color='#2f6fac', label='LLSG Contrast Curve')
	plt.legend(handles=[LLSG_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(LLSG_Distance, LLSG_Sensitivity_Gauss, linewidth =2, color='#2f6fac')
	
	#Saves the figure and then shows the plot. 
	savefig('new_LLSG_curves_{star_name}.png'.format(star_name=name_input))

	plt.show(Curve)
	
	#PCA contrast curve. 
	
	fig = plt.figure(figsize=(8,4))		#Initialises the plot.
	ax1 = fig.add_subplot(111)			#Creates the axis.
	
	#Creates and positions the legend.
	handles, labels = ax1.get_legend_handles_labels()
	line_up, = plt.plot([1,2,3], label='Line 2')
	line_down, = plt.plot([3,2,1], label='Line 1')
	plt.legend(handles=[line_up, line_down])
	ax1.legend(handles, labels)
	
	#Formats the colour and text for the legends, colour represented by hexadecimal.
	PCA_Legend = mpatches.Patch(color='#fcb141', label='PCA Contrast Curve')
	plt.legend(handles=[PCA_Legend])
	
	plt.xlabel('Angular separation [arcsec]')	#X Label.
	plt.ylabel(str(sigma)+' sigma contrast')	#Y Label.
	
	#Sets the y axis scale and the range of limits. 
	ax1.set_yscale('log')
	plt.ylim((0.000001,0.01))
	
	#Creates the variables that will be plotted.
	Curve = 0	#Initialises the Curve variable.
	Curve = plt.plot(PCA_Distance, PCA_Sensitivity_Gauss, linewidth =2, color='#fcb141')
	
	#Saves the figure and then shows the plot. 
	savefig('new_PCA_curve_{star_name}.png'.format(star_name=name_input))
	plt.show(Curve)


"""
Asks for positions of candidates/reference stars??
"""
#Input the position(x,y) and flux as array
def candidate(number):
	newcent = [512,512]
	ox = [newcent[0]]
	oy = [newcent[1]]
	cent_x = []
	cent_y = []
	position =[None]*num
	cent_x = ox*num
	cent_y = oy*num
	x = list()
	y = list()
	#f = list()
	m = None
	n = None 
	#l = None

	cent_x = ox*num
	cent_y = oy*num
	for i in range(int(num)):
		while True:
			m = input("Input the x co-ordinate of the reference/candidate position\n(Integer number)\n")

			while isinstance(m,int) == False:
				Loop_count = None
				print( "\nAt least one of these co-ordinates is not an integer.\nPlease enter an integer number\n" )
				m = input("Input the x co-ordinate of the reference/candidate position\n(Integer number)\n")

			if isinstance(m,int) == True: 
				break				

		while True:
			n = input("Input the y co-ordinate of the reference/candidate position\n(Integer number)\n")

			while isinstance(n,int) == False:
				print( "Please enter an integer number\n" )
				n = input("Input the y co-ordinate of the reference/candidate position\n(Integer number)\n")

			if isinstance(n,int) == True: 
				break	

		x.append(int(m))
		y.append(int(n))
		
	for i in range(int(num)):
		position [i] = (x[i],y[i])

	return position


"""
Converts from Cartesian to Polar coordinates
"""
def cart2pol (x,y,flux):

	rho = np.sqrt( x**2 + y**2)
	phi = np.degrees( np.arctan2(y,x) )
	phi = np.mod(phi,360)
	f_p = flux

	return rho,phi,f_p
	
	



"""
Function takes the cube and the PSF (plus its centre)
and fits a 2d Gaussian into it.
Using that, it recenters the cube and plots the result.
Returns the re-centered cube, and the x/y shifts
"""

def Gaussian_2d_Fit( psf, cube_orig, centre):

	gauss = vip.var.fit_2dgaussian(
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
		model='gauss', nproc = 1, subi_size = fwhm_int, 
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

print("\n\n\n ================= Project Script =================")
print("\n Version 11_py3 \n\n")



	
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
	------ 4 - RECENTER LOOP
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
		
		star_xy[1]=input("Input the x co-ordinate of the central position of the star (integer number): ")
		star_xy[1] = int(star_xy[1])
		
		
		#Opens the cube (using VIP) with cube_orig as HDU:0 and calls it cube_orig, 
		#and the parallactic angles from a text file.
		#Uses VIP to also open the point spread function previously loaded in.
		
		cube_orig = vip.fits.open_fits(cube)
		angs = np.loadtxt( '{name}_angles.txt'.format(name = name_input) )
		psf = vip.fits.open_fits(initialpsf, n=0, header=False, ignore_missing_end=True, verbose=True)
	
		print( "Fitting a 2D gaussian to centre the images..." )
		cube1, shy1, shx1, fwhm = Gaussian_2d_Fit( psf, cube_orig, star_xy)
		
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

	im1 = vip.preproc.cosmetics.frame_crop(cube_orig[0], 1000, verbose=False)
	im2 = vip.preproc.cosmetics.frame_crop(cube[0], 1000, verbose=False)
	
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

reference_xy[1] = input("Input the x co-ordinate of the reference/candidate position (Integer number): ")
reference_xy[1] = int(reference_xy[1])
	


"""
------ 6 - PCA
"""
print( 'Will now reduce the data using full frame PCA' )

#Reduce the image using full frame PCA then plots the image..
# Optimal number of PCs are calculated

ffpca_output = vip.pca.pca_fullfr.pca(
	cube, angs, fwhm = fwhm,
	source_xy = (reference_xy[0], reference_xy[1]),
	mask_center_px = None, ncomp = (1, len(cube) ,1),
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
	
	Rank = input("What rank of llsg would you like to do? (default =6)\n")
	Rank = int(Rank)
	
	#Reduces the image using the LLSG algorithm then plots the image.
	print(fwhm)
	fr_llsg = vip.llsg.llsg(cube, angs, fwhm=fwhm, rank = Rank, thresh = 2, max_iter = 20, random_seed = 10)
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
starphot = StarphotCalculation(Number_images, psf_xy, starphot)
starphot_array = np.zeros(500)

for i in range (0,500):
	starphot_array [i] = starphot

np.savetxt('{name}_starphot'.format(name=name_input), starphot_array)


"""
------ 9 - INJECT FAKE PLANETS
"""

"""
number_planets = 0.1 #Initialises the variable as an odd number to use the integer check loop.

#Loop checks that the number of planets is in fact an integer.
while isinstance(number_planets, int) == False or number_planets > 6:

	if Loop_count == 0:
		number_planets = input("Enter the integer number (1-6) of artificial planets to inject (recomended 5).\n")
	if Loop_count > 0:
		print( "\nThis is not an allowed value.\nPlease try again\n" )
		number_planets = input("Enter the number of artificial planets to inject.\n")

	Loop_count = Loop_count + 1 

#Sets loop count to 0 to be used again.
if isinstance(number_planets, int) == True:
	Loop_count = 0

#Loop for if artificial planets are used, loop calculates where to place them,
#injects them and then removes them again.
if number_planets > 0:
	rad_dists = np.empty(number_planets)

	#Calculates the orbital radius to place the injected planets at.	
	for i in range(number_planets):
		rad_dists[i] = seperation * (i+1)
	
	print( rad_dists )

	#Loop to check that the number of branches entered is an integer number.
	n_branches=0.1 #Initialises the variable to non integer so loop will run
	while isinstance(n_branches, int) == False:

		if Loop_count == 0:
			n_branches=input("How many branches of synthetic planets would you like to use (recomended 3)?\n")
		if Loop_count > 0:
			print( "\nThis is not an integer.\nPlease enter an integer number\n" )
			n_branches = input("How many branches of synthetic planets would you like to use?\n")

		Loop_count = Loop_count+1 

	#Sets Loop_count to 0 for use again.
	if isinstance(n_branches, int) == True:
		Loop_count = 0

	#Calculates the locations of the branches with their respective injected planets.
	if n_branches > 0:
		theta = 360 / n_branches	#The angle to separate each branch by.

	#Checks if the user would like to input their own 
	#brightness value for the injected planets.
	brightness_check=0
	Loop_count = 1 #allows error check of choice of flux input
	print( "1. Enter planet brightness manually (eg. flux of candidate)\n" )
	print( "2. Use background brightness multiplied by user value (recomended 3)\n" )

	while Loop_count:
		brightness_check = input("Enter choice (1,2):\t")
		if brightness_check == 1 or brightness_check == 2:
			Loop_count = 0
		else:
			print( "\nInvalid entry\n" )

	#Loop for the user to use their own value of flux.
	if brightness_check == 1:
		flvl = input("Please input the planets brightness here\n")

	#Loop for when the user does not want to enter a flux value and ask for multiple value of of background flux.
	if brightness_check == 2:
		Loop_count = 1 #allows loop to check if correct value entered ie not zero, loop breaks when Loop_count=0
		flux_multiple = 3
		while Loop_count:
			flux_multiple = input("Enter multiple value for background flux (3 recomended)\n")
			if flux_multiple < 0.00001:
				print( "\nInvalid Entry\n" )
			else:
				Loop_count = 0	
						
		averageflux = backgroundflux(cube, averageflux)
		flvl = float(averageflux)*float(flux_multiple)

#Creates a normalised PSF and a data cube with the injected planets.
psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=100, threshold=None, mask_core=None)
cubefc = vip.metrics.fakecomp.cube_inject_companions(cube, psf, angs, flevel=flvl, plsc=pxscale_keck,rad_dists=rad_dists,theta=theta, n_branches=n_branches) 

#Initialises a bunch of variables to do with calculating the location
#of the injected planets.

average = np.zeros(2)
sum = np.zeros(2)
newcent = np.zeros(2)
synplanetlocation = np.zeros(2)

#From the shift of frames in pixels, find the center of the image
Nx = np.prod(shx1.shape)
for i in range(0,Nx):
	sum[0] = sum[0] + shx1[i]
	
average[0] = sum[0]/Nx	

Ny = np.prod(shy1.shape)
for i in range(0,Ny):
	sum[1] = sum[1] + shy1[i]
	
average[1] = sum[1]/Ny	

newcent[0] = star_xy[0] + average[0] 
newcent[1] = star_xy[1] + average[1]


synplanetlocation[0] = newcent[0] + (seperation * sin(-theta))
synplanetlocation[1] = newcent[1] + (seperation * cos(-theta))
synplanetlocation[0] = round(synplanetlocation[0])
	#There will be a noticeable error if the shift between images is large
synplanetlocation[1] = round(synplanetlocation[1])
	#These are only rough guesses, tends to be out by a couple of pixels

print( "synplanetlocation = ", synplanetlocation)

# Optimises the number of principle components with the new cubefc.
opt_pcs = vip.pca.pca_optimize_snr(
	cubefc, angs, fwhm=fwhm, 
	source_xy=(synplanetlocation[0], synplanetlocation[1]),
	mask_center_px=None, fmerit='mean', range_pcs=(1,20)
	)

#Plots the data cube with the injected planets
fr_pca3 = vip.pca.pca(cubefc, angs, ncomp=opt_pcs)

vip.fits.write_fits('new_Planet_injectionPCA.fits', fr_pca3, verbose=True) 

_ = vip.metrics.frame_analysis.frame_quick_report(fr_pca3, fwhm=fwhm, source_xy=(synplanetlocation[0],synplanetlocation[1]))
	#give nan in here (?)
"""
	
"""
Some code commented here (planet subtraction)
-----------------------------------------------

source_xy=[(reference_xy[0], reference_xy[1])]

print( "source_xy = ", source_xy )

r_0, theta_0, f_0=firstguess(cubefc, angs, psf, ncomp=10, plsc=pxscale_keck,planets_xy_coord=source_xy,
fwhm=fwhm,annulus_width=3, aperture_radius=3,f_range=np.linspace(100,2000,20),  simplex=True, display=True, verbose=True)

# Planet subtraction
plpar_bpicb=[(r_0, theta_0, f_0)]
cube_emp=cube_planet_free(plpar_bpicb, cubefc, angs, psf, pxscale_keck)
fr_pca_emp=vip.pca.pca(cube_emp, angs, ncomp=opt_pcs, verbose=False)
pt_sub = 0.1
-----------------------------------------------

"""




"""
------ 10 - CONTRAST CURVES
"""

psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=101, threshold=None, mask_core=None)

PCA_contrast_curve = vip.metrics.contrcurve.contrast_curve(
	cube, angs, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=1, algo=vip.pca.pca, 
	ncomp=40, debug=True,save_plot='PCA'
	)

LLSG_contrast_curve =  vip.metrics.contrcurve.contrast_curve(
	cube, angs, psf, fwhm, pxscale_keck, starphot,
	sigma=sigma, nbranch=1, algo=vip.llsg.llsg, 
	debug=True,save_plot='LLSG'
	)

Contrastcurvedata(PCA_contrast_curve, LLSG_contrast_curve, pxscale_keck, sigma,name_input)



"""
------ 11 - STIM MAPS
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
------ 12 - ANDROMEDA MAP
"""




print( "===========================	End of the program	===========================" )










# Some space here