"""
============================	LIBRARIES	=================================
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import vip_hci as vip
import hciplot
import pandas as pd
from astropy.io import fits
import warnings
warnings.simplefilter("ignore")

#Additional functions from other .py files 
from ProjectScript_ultils import *	# import all functions
from ProjectScript_preprocessing import *	# import all functions

"""
============================	GLOBAL VARIABLES		=================================
"""
psf_xy = [95,111]			#Centre of psf.
reference_xy = [0.1, 0.1]	#Reference bright point defined as a decimal so it loops through the decimal check.
star_xy = [0.1, 0.1] 		#Centre of the star co-ordinates, defined as a decimal so it loops through the decimal check.
averageflux = 0				#Initialises the background flux.				
flvlup = 15000				#Just a range of flux's for the first guess algorithm .
flvldown = 1000				#to loop over.

#Initialises the loop counter, loop counter is used to prevent error warning on first trial
#when looping through integer check functions.
Loop_count = 0				
pxscale_keck = 0.00953 		#Pixel scale for keck 2 (can be found on the Keck Website)
seperation = 60 			#Radial seperation between synthetic planets in terms of pixels.
starphot = 0				#Initialises the Starphot parametre.
sigma = 5					#Sets the value of sigma for the contrast curves.
sub_count = 0
mask_center = 10			#Mask the central region
path_psf = 'psf.fits'

"""
============================	FUNCTION DEFINITIONS		=================================

	This script only contains functions for postprocessing, e.g. PCA (Full-frame, Annular), LLSG. 
	Detection module is avaliable from another script

"""

def FullFramePCA(cube, angs, fwhm, name_input):	
"""
Adapted from vip-hci. Returns FFPCA frame and optimal number of principal components for input image cube. 

Parameter: 
	cube : numpy ndarray, 3d
		Input ADI cube.
	angs : numpy ndarray, 1d
		Corresponding parallactic angle for each frame. Extracted from header of each .fits file
	fwhm : float
		Known size of the FHWM in pixels to be used.
	name_input : character string
		Name of the star (obtained from the first input)

Return:
	If full_output == True and source_xy != None, returns 
		0 -	Residuals cube : numpy.ndarray, 3d,(N,1024,1024) (N= number of images)
		1 -	Residual image : numpy.ndarray, 2d, (1024,1024)
		2 - Pandas dataframe with PCs, S/Ns and fluxes
			<class 'pandas.core.frame.DataFrame'> 10 (columns = PCs, S/Ns, fluxes)	
"""

	print( "\n--- FULL FRAME PCA")
	print( "\nUse DS9 to obtain the reference coordinates" )
	print( "\nA 'reference' is simply any kind of stable bright spot across the images, like a static speckle or a candidate.")

	reference_xy[0] = input("Input the x coordinate of the reference: ")
	reference_xy[0] = int(reference_xy[0])

	reference_xy[1] = input("Input the y coordinate of the reference: ")
	reference_xy[1] = int(reference_xy[1])

	print( '\nApplying Full Frame PCA...' )

	ffpca_output = vip.pca.pca_fullfr.pca(
				cube, angs, fwhm = fwhm,
				source_xy = (reference_xy[0], reference_xy[1]),
				mask_center_px = mask_center, ncomp = (1, len(cube) ,1),
				full_output=True)

	# Save FFPCA output
	residual_cube = ffpca_output[0]
	fr_pca1 = ffpca_output[1]

	# Getting Optimal Number of PCs
	# Get pandas dataframe table output into a numpy array
	pca_data = ffpca_output[2].rename_axis('ID').values

	# Extract and save S/N ratio and PCs columns in arrays
	snr_data = []
	pcs_data = []

	for i in range(0, len(pca_data)):
		pcs_data.append(pca_data[i][0])
		snr_data.append(pca_data[i][1])

	# Get the index of the maximum value of the S/N column, 
	# and retrieve that same value in the PCs column
	# This will be the optimal number of principal components
	snr_max = np.argmax(snr_data, axis=0)
	optimal_pcs = int(pcs_data[snr_max])

	print("Optimal PCs", optimal_pcs)

	# Loop asks user if they would like to save the image.
	print( 'Would you like to save the Full Frame PCA reduced images?' )
	questionsave = Checkfunction()

	if questionsave == 0:
		vip.fits.write_fits('new_FFPCA_{name}.fits'.format(name=name_input), fr_pca1, verbose=True) 

	return residual_cube, fr_pca1, optimal_pcs


def AnnularPCA(cube, angs, fwhm, name_input):
"""
Adpated from vip_hci. Returns annular pcs frame, and residuals cube for STIM map. 

Parameter: 
	cube : numpy ndarray, 3d
		Input ADI cube.
	angs : numpy ndarray, 1d
		Corresponding parallactic angle for each frame. Extracted from header of each .fits file
	fwhm : float
		Known size of the FHWM in pixels to be used.
	name_input : character string
		Name of the star (obtained from the first input)

Return: 
	annpca_frame: numpy ndarray, 2d
		Residual image

	residuals_stim_input: numpy ndarray, 3d
		Residual image 

"""

	print( "\n --- ANNULAR PCA")

	print("Performing Annular PCA...")

	ann_pca_output = vip.pca.pca_local.pca_annular(
							cube, angs, fwhm=fwhm,asize=fwhm/2,radius_int=mask_center,
							full_output = True, verbose = True
							)

	# Outputs three objects in this order:
	#	cube_out: cube of residuals
	#	cub_der: derotated cube of residuals
	#	frame: Annular PCA frame
	
	vip.fits.write_fits('AnnPCA_{name}.fits'.format(name=name_input), ann_pca_output[2], verbose=True)
	vip.fits.write_fits('AnnPCA_residuals_{name}.fits'.format(name=name_input), ann_pca_output[1], verbose=True) 

	annpca_frame = ann_pca_output[2]
	residuals_stim_input = ann_pca_output[1]

	return annpca_frame, residuals_stim_input 
	#change residuals_stim_input to annpca_residual_cube

def LLSG(cube, angs, fwhm, name_input):
"""
-----------------	LLSG	-----------------

Returns LLSG frame, and residuals cube for STIM map. 

Parameter: 
	cube : numpy ndarray, 3d
		Input ADI cube.
	angs : numpy ndarray, 1d
		Corresponding parallactic angle for each frame. Extracted from header of each .fits file
	fwhm : float
		Known size of the FHWM in pixels to be used.
	name_input : character string
		Name of the star (obtained from the first input)

Return:	
	llsg_frame : numpy ndarray, 2d
		Residual frame

	llsg_residual_sparse : numpy ndarray, 3d
		Residual Cube

	Rank : integer 
		Rank of 



	

"""

	print( "\n   --- LLSG")
	Rank = input("What rank of llsg would you like to do? (default =6)\n")
	Rank = int(Rank)

	#Reduces the image using the LLSG algorithm then plots the image.
	print(fwhm)

	llsg_output = vip.llsg.llsg(cube, angs, rank = Rank, fwhm=fwhm, thresh = 2, max_iter = 20, random_seed = 10,full_output = True)

	"""
	llsg_output = vip.llsg.llsg(cube, angs, rank = None, fwhm=fwhm, 
							max_iter = len(cube)-1, full_output = True, auto_rank_mode='noise',
							low_rank_mode='svd', residuals_tol=1e-1, cevr=0.9, thresh_mode='soft', nproc=1, n_segments=4, azimuth_overlap=None, radius_int=mask_center,
							random_seed=None, imlib='opencv', interpolation='lanczos4',
							high_pass=5, collapse='median', verbose=True,
							debug=True)
	"""

	# LLSG full output is:
	#0	list_l_array_der
	#1	list_s_array_der
	#2	list_g_array_der
    #3	frame_l
    #4	frame_s
    #5	frame_g
    # Each is the frame and residual cube for the three noise types:
    # l (low-rank), s (sparse), g (gaussian)
    # Companions are in the Sparse noise.
    # To get STIM map of LLSG, input 'list_s_array_der' into STIM Map

	llsg_frame = llsg_output[4]
	llsg_residual_sparse = np.asarray( llsg_output[1][0] )
		#list_s_array_der is a list of one member,
		# this one member is another list of the 24 residual images
		# So get that one element, and make it an array.

	hciplot.plot_frames(llsg_frame, label ='LLSG reduced image of {name}'.format(name=name_input))
	
	
	#Loop asks user if the would like to save the image.
	print( 'Would you like to save the LLSG reduced image?\n' )
	questionsave = Checkfunction()

	if questionsave == 0:
		vip.fits.write_fits('new_LLSG_{name}.fits'.format(name=name_input), llsg_frame, verbose=True)


	return llsg_frame, llsg_residual_sparse, Rank

"""
-----------------	STIM MAP	-----------------
"""
def StimMap(residuals_cube, name_input, origin_algo):

	print("Computing STIM map using {algo} output...".format(algo = origin_algo))
	stim_map = vip.metrics.compute_stim_map(residuals_cube)


	hciplot.plot_frames( 
			stim_map, 
			label = 'STIM map of {name} using {algo}'.format(name=name_input, algo = origin_algo), 
			grid = False, 
			size_factor = 5
			)

	#Loop asks user if they would like to save the image.
	print( 'Would you like to save the image?' )
	questionsave = Checkfunction()

	if questionsave == 0:	
		vip.fits.write_fits('STIM_{algo}_{name}.fits'.format(algo=origin_algo, name=name_input), stim_map, verbose=True)


"""
-----------------	CONTRAST CURVES	-----------------
"""
def ContrastCurves(cube, angs, psf, fwhm, starphot, optimal_pcs, name_input):


	print(" Performing Starphot calculation...")
	#Starphot is the  total aperture flux from the star during
	#the integration time (if no coronograph).
	starphot = StarphotCalculation(name_input, psf_xy, starphot, fwhm)

	print(" Building contrast curves...")
	psf = vip.var.shapes.get_square(psf, 101, y=psf_xy[1], x=psf_xy[0], force=True, verbose=False)

	AnnPCA_contrast_curve = vip.metrics.contrcurve.contrast_curve(cube, angs, psf, fwhm, pxscale_keck, starphot, fc_snr = 10, cube_ref=None, 
		sigma=sigma, algo=vip.pca.pca_local.pca_annular, full_output = True, save_plot='RDI_curve', 
        radius_int=mask_center, asize=fwhm/2, ncomp='auto',nproc=1)
#
	PCA_contrast_curve = vip.metrics.contrcurve.contrast_curve(
		cube, angs, psf, fwhm, pxscale_keck, starphot,fc_snr = 10, 
		sigma=sigma, nbranch=1, algo=vip.pca.pca, 
		ncomp=optimal_pcs, debug=True,save_plot='PCA_curve', mask_center_px = mask_center)


	#for i in range (11,40):
	#	LLSG_contrast_curve =  vip.metrics.contrcurve.contrast_curve(
	#		cube, angs, psf, fwhm, pxscale_keck, starphot,fc_snr = 10, 
	#		sigma=sigma, nbranch=1, algo=vip.llsg.llsg, 
	#		debug=True, save_plot='LLSG_curve', rank = i, radius_int=mask_center)
#
	#	outputLLSG = name_input + "_LLSG_output" + "_{rank}".format(rank = i)
	#	LLSG_contrast_curve.to_csv('{file}.txt'.format(file=outputLLSG), sep='\t', index=False)


	output_name1 = name_input + "_fc_all.fits"
	output_name2 = name_input + "_nofc.fits"
	##Saving all outputs
	vip.fits.write_fits(output_name1, AnnPCA_contrast_curve[1], verbose=True) 
	vip.fits.write_fits(output_name2, AnnPCA_contrast_curve[2], verbose=True)

	#Saves the PCA,ADI and LLSG curve outputs
	#Converting panda frame into .txt file
	outputAnnPCA = name_input + "_AnnPCA_output"
	outputPCA = name_input + "_PCA_output"
	#outputLLSG = name_input + "_LLSG_output"

	AnnPCA_contrast_curve[0].to_csv('{file}.txt'.format(file=outputAnnPCA), sep='\t', index=False)
	PCA_contrast_curve.to_csv('{file}.txt'.format(file=outputPCA), sep='\t', index=False)
	#LLSG_contrast_curve.to_csv('{file}.txt'.format(file=outputLLSG), sep='\t', index=False)












"""
-----------------	ANDROMEDA	-----------------
"""
def Andromeda(cube, angs, psf, fwhm):
	
	print(" Running ANDROMEDA algorithm...")
	print(" (may take an hour or more to complete)")
	
	# Normalises PSF (Done above too)
	psf = vip.metrics.fakecomp.normalize_psf(psf, fwhm, size=100,
										threshold=None, 
										mask_core=None,
										force_odd = False)

	# Andromeda requires even-sized images, hence "force_odd = False" is used.


	# In NIRC2, the Kp filter has central wavelength 2.124 um
	#			and Bandpass width 0.351 um
	# Source: https://www2.keck.hawaii.edu/inst/nirc2/filters.html

	print("Cube has shape", cube.shape)
	print("PSF has shape", psf.shape)

	survey_wavelength = 2.124
	shannon_wavelength = 0.351

	oversampling_factor = survey_wavelength/shannon_wavelength

	# Ran code on HIP544 with default settings and iwa=1. Outputs:
	"""
	WARNING: 18 frame(s) cannot be used because it wasn't possible to find any other frame to couple with them. 
	Their indices are: [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
	For all frames to be used in this annulus, the minimum separation must be set at most to 0.3631489127511684 *lambda/D 
	(corresponding to 4.395032995347474 pixels).
	"""
	# The set min_sep to 0.36


	andromeda_output = vip.andromeda.andromeda(
						cube,
						oversampling_factor,
						angs,
						psf,
						iwa = 1,
						min_sep = 0.36,
						verbose = True
						)

	print("Andromeda output length", len(andromeda_output))

	contrast_andr = andromeda_output[0]
	snr_map_andr = andromeda_output[1]
	snr_norm_map_andr = andromeda_output[2]
	stdcontrast_map_andr = andromeda_output[3]
	stdcontrast_norm_andr = andromeda_output[4]
	likelihood_andr = andromeda_output[5]
	ext_radius_andr = andromeda_output[6]

	vip.fits.write_fits('ANDR_contrast_{name}.fits'.format(name=name_input), contrast_andr, verbose=True)
	vip.fits.write_fits('ANDR_likelihood_{name}.fits'.format(name=name_input), likelihood_andr, verbose=True)

	print("ext_radius = ", ext_radius_andr)

		
"""
******************************************************************
				Main Body of the Project Script
******************************************************************
"""

# ****************************************************************
#						PRE PROCESSING
# ****************************************************************


print("\n\n\n ================= Project Script =================")
print("\n Version 13 \n\n")
	
"""
------ 1 - STAR NAME AND PRE-PROCESSING LOOP
"""
	
#name_input = input('Name of the star (example = HIP_544):   ')

name_input = 'BD+45598'

if (name_input[0] == 'B'):
	dirname = name_input
	nameTemp = name_input.split("+")
	filename = nameTemp[0] + '_' +nameTemp[1]
else: 
	nameTemp = name_input.split("_")
	dirname = nameTemp[0] + nameTemp[1]

print("\n\nName of the directory: ",dirname)
os.chdir('{name}/'.format(name=dirname))
os.system('ls')
#epoch = input('Epoch of the star from listing (example = 2012aug28):  ')
epoch = '2014oct03'

os.chdir('{date}/'.format(date=epoch))
Centring_loop = 0 # Initialise loop to allow for re centring

print("\nImport centeredcube?\n")
#cube_input = Checkfunction()
cube_input = 0
name_input = filename

if cube_input == 0:
	#Import the centred cube
	path_centcube = 'centeredcube_'+ name_input + '.fits'
	hdulist = fits.open(path_centcube, ignore_missing_end=True,verbose=False)
	cube_orig = hdulist[0].data
	cube = vip.preproc.cosmetics.cube_crop_frames(cube_orig, 500, xy=(511,511), force=True, verbose=True, full_output=False)
	psf = vip.fits.open_fits(path_psf, n=0, header=False, ignore_missing_end=True, verbose=False)
	fwhm = Calculate_fwhm(psf,psf_xy)
	angs = np.loadtxt('{name}_angles.txt'.format(name=name_input))

else: 
	cube, fwhm, angs, psf = CenteringCube()

# ****************************************************************
#						POST PROCESSING
# ****************************************************************

# Full Frame PCA
ffpca_residuals, ffpca_frame, optimal_pcs = FullFramePCA(cube, angs, fwhm, name_input)

# Annular PCA
#annpca_frame, pca_residuals = AnnularPCA(cube, angs, fwhm, name_input)

# LLSG
#llsg_frame, llsg_residuals, Rank = LLSG(cube, angs, fwhm, name_input)

# STIM Map
#StimMap(ffpca_residuals, name_input, "FFPCA")
#StimMap(pca_residuals, name_input, "PCA")
#StimMap(llsg_residuals, name_input, "LLSG")


# Contrast Curves
ContrastCurves(cube, angs, psf, fwhm, starphot, optimal_pcs, name_input)
#Contrastcurvedata(pxscale_keck, sigma, name_input)

# Andromeda
# Andromeda(cube, angs, psf, fwhm)
# good, bad = vip.stats.clip_sigma.clip_array(PCA, lower_sigma = 0, upper_sigma = 1000, out_good=False, neighbor=False,num_neighbor=None, mad=False)
# final = vip.stats.clip_sigma.sigma_filter(PCA, bad, neighbor_box=3, min_neighbors=3, verbose=False)

print( "===========================	End of the program	===========================" )

