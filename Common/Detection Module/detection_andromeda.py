

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import hciplot
import math
import pandas as pd
import vip_hci as vip
import os


from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm

"""
	Python implementation of IDL - Andromeda Detection Module
	Original Author: F. Cantalloube
	Implementation by: Jorge Fernandez (jorgefz.fernandez@gmail.com)

	Implemented from version:
		Revision: 1.5, Date: 2018/06/25 
"""


"""
	Andromeda Detection Module
	--------------------------

	Procedure to automatically detect companions in the SNR-map provided
    by ANDROMEDA and derive its flux in the contrast-map provided by ANDROMEDA.
    The detections are sorted out and displayed, directly on screen 
    (images and info) and both displayed images and information about the 
    detections (in a data file) are recorded.
    This procedure optionaly computes and displays the detection limit
    of the image set when processed by ANDROMEDA.
    A certain number of information are then output to inform about the
    whole process, see paper Cantalloube et al. (2015) for more info.

	
	ARGUMENTS :

---Required input to detect and sort out likely planetary signals:

	snr_map				: 2D float array
						SNR-map provided by ANDROMEDA
						-> For detection and subpx astrometry.

	snr_norm			: 2D float array, optional.
						Normalised SNR map provided by ANDROMEDA

	snr_stddev			: 2D float array, optional.
						Standard Deviation SNR map.

	flux_map 			:  2D float array  
						FLUX-map provided by ANDROMEDA.
						-> For photometry.

	flux_stddev			: 2D float array
                        Map of the standard deviation  
                        of the estimated flux provided by ANDROMEDA.
                        -> For error on contrast estimation
                        -> And for detection limit computation.

	flux_norm_stddev	:2D float array
						Map of normalised standard deviation of flux.

	threshold			: integer (default = None), optional.
						threshold chosen so that detections are signal.
						found in the SNR-map that are above this threshold.

	size_subwindow		: even integer [pixels]
						size of the 'subimages' fully containing the signal
						and in which the 2D-Gaussian fit are performed.

	pixscale			: float [arcsec/pixelSize]
						equivalent size of one pixel in arcsec 
						(given by the camera used for imagery).

	oversampling		: float [pixels]
						oversampling of the images given by the camera.

	max_lobe_dist		: integer [pixels], optional
						maximum distance between the main lobe and the
						tertiary (positive) lobe in the pattern 
						of a planet signature made in the SNR-map. 
						If the main lobe has a very high SNR, it 
						may happen that the tertiary lobe is above 
						the threshold and thus detected as a companion
						whilst being an artefact.

	owa					: float [lambda/D], optional
						External limit of the data in the images (radius)

	iwa					: integer [lambda/D], optional
						Internal limit of the data in the images (radius),
						from which planetary signals are sought.

	obs_flag			: integer, optional
						Defines constraints for Gaussian Fit
						0 - NaCo
						1 - NIRC2
						2 - STIM map

	save_plots			: bool, optional
						If true, generates 3d plots of signals and their gaussian fits
						on 'detections' folder.

---Required input to compute the contrast in astrophysical terms:

	flux_map   			: 2D float array. Flux-map provided by ANDROMEDA.
	psf 			 	: 2D float array. PSF used to launch ANDROMEDA.
	img_int_time		: float. Integration time of the images.
	psf_int_time		: float. Integration time of the PSF.
	tnd					: (?). Transmission of the neutral density, 
						if one used to image the PSF that acts as a reference.

---Required input to save information about the detection process:

	detection_files		: string array 
						Name of the files that are saved,
						containing the properties of all the detected signals.
						(.dat, .fits, .eps)

	nb_detections		: (output) - integer
						Number of probable detections made.

---Others:

	verbose				: bool (default = True), optional.
						Displays coments while compiling the detection.

	naco				: bool (default = True), optional.
						True if analysing NaCo data.


	IDL Functions Required:
		im_contok.pro
		mpfit2dpeak.pro (and all subfunctions needed)
		gauss_calc.pro (cantallf)
		saveimage.pro 
		detection_limit.pro (cantallf)
		tvwin.pro (ONERA);   
	

	"""




# Main class that stores maps
class Maps:
	def __init__(self, snr_map, snr_norm, snr_stddev,
			flux_map, flux_stddev, flux_norm_stddev):
		# SNR maps
		self.snr_map = snr_map
		self.snr_norm = snr_norm
		self.snr_stddev = snr_stddev
		# Flux maps
		self.flux_map = flux_map
		self.flux_stddev = flux_stddev
		self.flux_norm_stddev = flux_norm_stddev
		#Flags & other variables
		self.flux_flag = True
		self.snrnorm_flag = True

	def check_input(self, verbose):
		# Checking SNR maps
		if (not isinstance(self.snr_map, np.ndarray)):
			raise ValueError(" Error: snr map is not a numpy array.")
		elif(self.snr_map.ndim != 2):
			raise ValueError(" Error: snr_map must have 2 dimensions.")

		if not isinstance(self.snr_stddev, np.ndarray):
			if verbose:
				print(" No Standard Deviation SNR map provided. Will be generated instead.")
			self.snr_stddev = np.ones((self.snr_map.shape[0], self.snr_map.shape[0]), dtype = float)

		if not isinstance(self.snr_norm, np.ndarray):
			if verbose:
				print(" No normalised SNR map provided. Non-normalised one will be used instead.")
			self.snr_norm = np.copy(self.snr_map)
			self.snrnorm_flag = False

		# Checking Flux maps
		if not isinstance(self.flux_map, np.ndarray):
			if verbose:
				print(" No correct Flux map provided. No flux calculations will be performed.")
			self.flux_flag = False

		if (self.flux_flag == True) and (not isinstance(self.flux_stddev, np.ndarray)):
			if verbose:
				print(" No Standard Deviation Flux map provided. Will be generated instead.")
			self.flux_stddev = np.ones((self.snr_map.shape[0], self.snr_map.shape[0]), dtype = float)
		
		if (self.flux_flag == True) and (not isinstance(self.flux_norm_stddev, np.ndarray)):
			if verbose:
				print(" No Normalised Standard Deviation Flux map provided. Will be generated instead.")
			self.flux_stddev = np.ones((self.snr_map.shape[0], self.snr_map.shape[0]), dtype = float)


# Man class that stores important variables other than maps
class Data:
	def __init__(self, threshold, pixscale, oversampling, owa, iwa, 
					size_subwindow, max_lobe_dist, obs_flag, size_snr):
		self.threshold = threshold
		self.pixscale = pixscale
		self.oversampling = oversampling
		self.owa = owa
		self.iwa = iwa
		# Size of subwindow from which to cut-out signals
		self.sz = size_subwindow
		# Size of input science images in pixels
		self.img_size = size_snr
		# Min distance between two signals to be considered separate
		self.max_lobe_dist = max_lobe_dist
		# Observation flag. Fit constraints for NaCo (0), NIRC2 (1), or STIM (2)
		self.obs_flag = obs_flag
		# Size of element of resolution in lambda/D
		self.res = 1.0
		# size of the subwindow in lambda/D
		self.sz_res = 1.0

	def check_input(self, verbose, flux_ns):
		#Checking pixscale
		if not isinstance(self.pixscale, (float,int)):
			raise ValueError(" No input pixscale or wrong type. Must be float or int greater than zero.")
		if (self.pixscale <= 0):
			raise ValueError(" Pixel scale must be greater than zero.")
		#Checking oversampling
		if not isinstance(self.oversampling, (float,int)):
			raise ValueError(" No input oversampling or wrong type. Must be float or int greater than zero.")
		if (self.oversampling <= 0):
			raise ValueError(" Oversampling factor must be greater than zero.")
		#Checking outer working angle (owa).
		# Looks for zero values from edge of flux_norm_stddev map towards the middle.
		# Defines owa as the distance from the first non-zero value to the centre
		if (not self.owa) or (isinstance(self.owa, (float,int)) and self.owa <= 0):
			"""
			mid = int((self.img_size-1)/2)
			slab = flux_ns[:,mid]
			mid2 = math.ceil(1 + slab.shape[0] / 2)
			sliced = slab[1:mid2]
			ind = 0
			for i,arr in enumerate(sliced):
				if (arr != 0):
					ind = i+1
					break
			"""
			self.owa = (self.img_size/2 - 1) / (2*self.oversampling)
			if verbose:
				print(" No input owa set or invalid value. Using %.3f lambda/D." % self.owa)
		#Checking inner working angle (iwa)
		if (not self.iwa) or (isinstance(self.iwa, (float,int)) and self.iwa <= 0):
			self.iwa = 1.0
			if verbose:
				print(" No iwa set or invalid value. Using %.2f. lambda/D." % self.iwa)
		# Checking subwindow size
		if (not self.sz) or (isinstance(self.sz, (float,int)) and self.sz <= 0):
			self.sz = math.floor(4*self.oversampling) + 2
			if verbose:
				print(" No subwindow size set or invalid value. Using %d." % self.sz)
		self.sz = 1.5*self.sz
		# Checking that subwindow size is odd
		if (self.sz % 2 == 0):
			self.sz += 1
			if verbose:
				print(" Sub-image size must be odd. Value set to %d." % self.sz)
		# Checking maximum main-tertiary lobe distance
		if (not self.max_lobe_dist) or (isinstance(self.max_lobe_dist, (float,int)) and self.max_lobe_dist <= 0):
			self.max_lobe_dist = math.floor(10*self.oversampling)
			if verbose:
				print(" No max lobe distance set or invalid value. Using %.2f." % self.max_lobe_dist)
		# Checking Threshold
		if (not self.threshold) or (isinstance(self.threshold, (float,int)) and self.threshold <= 0):
			if verbose:
				print(" No threshold input set or invalid value. Using %.2f sigma." % self.threshold)
			self.threshold = 5.0
		# Checking observation flag
		if not isinstance(self.obs_flag, int):
			self.obs_flag = 0
			if verbose:
				print(" Invalid observation flag. Using %d." % self.obs_flag)
		# Size of element of resolution in lambda/D
		self.res = 2*self.oversampling
		# size of the subwindow in lambda/D
		self.sz_res = float(self.sz)/self.res



class idl:

	def polaire2(rt, width = None, oc = 0.0,
				center_x = None, center_y = None,
				even_px = False, leq = False,
				verbose = True):

		#---------------------
		# 	Input Checking
		#---------------------

		# Checking RT:
		if not rt:
			print("Error: No input telescope radius!")
			return

		if not isinstance(rt, (float,int)):
			print("Telescope radius must be a float or int")
			return

		if (rt <= 0):
			print("Telescope radius must be greater than zero")
			return

		# Checking OC:
		if not oc:
			oc = 0.0

		if (oc < 0) or (oc > 1):
			print("Occultation must be in range [0,1]. Will be set to 0")
			oc = 0.0

		# Checking Width
		if not isinstance(width, (float,int)) or ( isinstance(width, (float,int)) and (width < 0) ):
			if even_px:
				width = 2 * idl.round(rt)
			else:
				width = 2 * int(rt) + 1

		# Definiton of center x,y
		if not center_x:
			center_x = float( (width-1.0)/2 )

		if not center_y:
			center_y = float( (width-1.0)/2 )

		# Verbose
		if verbose:
			print(" Width = ",width)
			print(" Center (x,y) = ", center_x, center_y)


		# Calculating Angles Rho and Phi
		width = int(width)

		x = np.arange(width*width, dtype=float)
		x = x.reshape((width,width))
		x = x % width
		# These 3 lines are equivalent to IDL:
		#	x = float( dindgen(widht,width) mod width )

		y = np.transpose(x)

		x = x - center_x
		y = y - center_y

		rho = np.sqrt(x**2 + y**2)/float(rt)
		phi = np.arctan2(y, x + (rho==0) )
		# Here (rho==0) returns a matrix of zeros (False) where the element of rho is not ==0
		# and True (1) is the element is ==0
		# Hence, x + (rho == 0), will add 1 to the element of X in which the corresponding
		# element in rho is ==0. 

		# Calculating Mask
		if leq:
			mask = ( (rho <= 1)*(rho >= oc) ).astype(int)
		else:
			mask = ( (rho < 1)*(rho >= oc) ).astype(int)

			# Get matrix (rho <= 1) (array of bools)
			# Get matrix (rho >= oc) (array of bools)
			# Multiply both. (Note: True*True = True; False*True = False)
			# Convert from bool to integer array using numpy.astype

		return rho, phi, mask

	# ECLAT Function from IDL
	# Shifts x,y values in 2D array by rows/2 and cols/2, respectively.
	def eclat(img):
		rows = img.shape[0]
		cols = img.shape[1]
		new_img = np.roll(img, int(rows/2), axis=0)
		new_img = np.roll(new_img, int(cols/2), axis=1)
		return new_img



	"""
	 Based on  distc.pro function from IDL
	 Original Author: L. Mugnier.
	 Creates square matrix of size 'n*m'
	 where the values are proportional to their distance to the origin.
	 INPUTS:
		n 	:(int), rows
		m 	:(int), columns, optional
		cx 	:(int), x-center of matrix, optional
		cy 	:(int), y-center of matrix, optional
	"""
	def distc(n, cols=None, cx=None, cy=None):

		n = int(n)
		m = cols

		# Without columns, make the matrix square
		if not m:
			m = n

		# Initialize the return array
		arr = np.ones((n,n))

		x = np.arange(n)
		
		# Without center coordinates:
		if (not cx) and (not cy):
			for i in range(n):
				if x[i] >= (n-x)[i]:
					x[i] = (n-x)[i]
			x = x**2

			i = 0
			while(i < m/2):
				y = np.sqrt(x + i**2)
				arr[i] = y
				if (i != 0):
					arr[m-i] = y
				i += 1


		# If at least one center coordinate is input:
		else:
			if not cx:
				cx = n/2.
			if not cy:
				cy = m/2.
			x = (x - cx)**2

			i = 0
			while (i < m):
				arr[i] = np.sqrt(x + (i - cy)**2)
				i += 1

		return arr



	# Implementation of WHERE function from IDL
	# Original Author: Ralf Farkas (github.com/r4lv)
	# Source: github.com/vortex-exoplanet/VIP/blob/master/vip_hci/andromeda/utils.py
	def where(array_expression):
		res = np.array([i for i, e in enumerate(array_expression.flatten()) if e])
		return res

	# Implementation of ROUND function from IDL
	# Original Author: Ralf Farkas (github.com/r4lv)
	# Source: github.com/vortex-exoplanet/VIP/blob/master/vip_hci/andromeda/utils.py
	def round(x):
		val = np.trunc(x + np.copysign(0.5, x))
		return int(val)


# This class performs a Gaussian 2D Fit on an input image
# and stores the return values
class Gauss:
	def __init__(self, image):

		gauss_out = vip.var.fit_2dgaussian( image,
					full_output = True, debug=False)
		"""
		 'gauss' is a pandas dataframe with columns:
			0-centroid_y
			1-centroid_x
			2-fwhm_y
			3-fwhm_x
			4-amplitude
			5-theta
		"""

		# Centroid
		self.cx = gauss_out.iloc[0,1]
		self.cy = gauss_out.iloc[0,0]
		
		#Full Width Half Maxima
		self.fwhm_x = gauss_out.iloc[0,3]
		self.fwhm_y = gauss_out.iloc[0,2]

		self.peak = gauss_out.iloc[0,4] 	# Peak Amplitude
		self.angle = gauss_out.iloc[0,5] 	# Angle of fit
		self.good_fit = True	# Bool for quality of fit


# This class stores the range of values of a 'good' Gaussian Fit 
# and checks whether an input fit satisfies them
class Gauss_params:

	def __init__(self, res, sz, obs_flag):

		self.sz = sz

		if (obs_flag != 0) or (obs_flag != 1) or (obs_flag != 2):
			obs_flag = 0

		# NACO (VLT)
		if obs_flag == 0:
			self.lsize = [ 0.95*res, 0.25*res ]		# Max-min range for large fwhm
			self.ssize = [ 0.65*res, 0.20*res ]		# Max-min range for small fwhm
			self.cmax = 0.75		# Max deviation from sz center
			#-----------------------
			self.fwhm_ratio = 2		# Max value of fwhm_x/fwhm_y or fwhm_y/fwhm_x
			self.fwhm_min = (0.2)*sz 	# Min FWHM, if gaussian is too thin	
			self.fwhm_max = 1.5*sz  	# Max FWHM, if gaussian is too wide

		# NIRC2 (Keck)
		elif obs_flag == 1:
			self.lsize = [ 0.60*res, 0.25*res ]
			self.ssize = [ 0.65*res, 0.15*res ]
			self.cmax = 1.50
			self.fwhm_ratio = 2
			self.fwhm_min = (0.2)*sz
			self.fwhm_max = 1.5*sz

		# STIM map
		elif obs_flag == 2:
			self.lsize = [ 0.60*res, 0.25*res ]
			self.ssize = [ 0.65*res, 0.15*res ]
			self.cmax = 1.50
			self.fwhm_ratio = 2
			self.fwhm_min = (0.2)*sz
			self.fwhm_max = 1.5*sz


	def check(self, gauss, verbose):

		# Check if fit center is outside window
		if (gauss.cx <= 0) or (gauss.cx >= self.sz):
			print(" 	->Gaussian fit: Error on x-position")
			gauss.good_fit = False

		if (gauss.cy <= 0) or (gauss.cy >= self.sz):
			print(" 	->Gaussian fit: Error on y-position")
			gauss.good_fit = False
		
		# Get large/small fwhm
		if (gauss.fwhm_x > gauss.fwhm_y):
			fwhm_large = gauss.fwhm_x
			fwhm_small = gauss.fwhm_y

		else:
			fwhm_large = gauss.fwhm_y
			fwhm_small = gauss.fwhm_x

		# Check if gaussian is too elongated
		if (fwhm_large/fwhm_small) > self.fwhm_ratio :
			print(" 	->Gaussian fit: Too elongated")
			gauss.good_fit = False

		# Check size of FWHMs (fit not too wide/thin)
		if (fwhm_large > self.fwhm_max):
			print("		->Gaussian fit: FWHM too wide")
			gauss.good_fit = False

		if (fwhm_small < self.fwhm_min):
			print("		->Gaussian fit: FWHM too thin")
			gauss.good_fit = False
		
		"""
		print(" ----> Checkign Gaussian Fit:")
		print(" Small FWHM", fwhm_small," must be within ", self.ssize[0],"and", self.ssize[1])
		print(" Large FWHM", fwhm_large," must be within", self.lsize[0],"and", self.lsize[1])
		print(" Fit Center", np.abs( (sz-1)/2. - gauss.cx)," and ", np.abs( (sz-1)/2. - gauss.cy), " > ", self.cmax,"?")

		if (fwhm_small >= self.ssize[0]) or (fwhm_small <= self.ssize[1]):
			gauss.good_fit = False

		if (fwhm_large >= self.lsize[0]) or (fwhm_large <= self.lsize[0]):
			gauss.good_fit = False

		if ( np.abs( (sz-1)/2. - gauss.cx) > self.cmax ):
			gauss.good_fit = False

		if ( np.abs( (sz-1)/2. - gauss.cy) > self.cmax ):
			gauss.good_fit = False
		"""

		# Print verbose message
		if (not gauss.good_fit) and (verbose):
			print(" Warning: Gaussian fit is not correct.")	





# Makes a 3D plot with a data image and its gaussian fit
def plot_3d_window_gaussian(sz, gauss, window, filename=None):

	# Checking filename input
	if filename is None:
		filename = "outfile"

	if not isinstance(filename, str):
		filename = str(filename)


	fig = plt.figure(figsize=plt.figaspect(0.5))

	# Declaring meshgrid
	x = np.linspace(0, sz, sz-1)
	y = np.linspace(0, sz, sz-1)
	x, y = np.meshgrid(x, y)

	# ----------------
	# 1st Plot: Window
	# ----------------
	ax = fig.add_subplot(1, 2, 1, projection='3d')

	surf = ax.plot_surface(x, y, window, cmap=cm.coolwarm)
	fig.colorbar(surf, shrink=0.5, aspect=10)

	# ----------------
	# 2nd Plot: Gaussian Fit
	# ----------------
	ax = fig.add_subplot(1, 2, 2, projection='3d')

	sigma_x = gauss.fwhm_x / ( 2*math.sqrt(2*math.log(2)) )
	sigma_y = gauss.fwhm_y / ( 2*math.sqrt(2*math.log(2)) )
	
	# Shift center
	sh = (sz-1)/2.

	fit = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((x-sh)**2/(2*sigma_x**2)
    	+ (y-sh)**2/(2*sigma_y**2))))

	surf2 = ax.plot_surface(x, y, fit, cmap='Blues')
	fig.colorbar(surf2)

	fig.savefig("detections/"+filename+".png")

	#plt.show()

	plt.close(fig)



# Implementation of POLAIRE2 function from IDL
# Original Author: L. Mugnier.
# Generates the polar coordinates RHO and PHI 
# and the MASK of intensity (pupil) of a telescope 
# with RT radius and central linear occultation rate OC

	"""
	PARAMETERS:
		rt 		: (float, int), positive int or float.
				Telescope radius (in pixels)
		leq 	: (bool), optional
				Is only used to define the Mask; 
				if this flag is set, the mask is defined by the 
				remote points <= RT (and> = OC * RT). 
				By default it is defined by the points to distance < RT
				(and >= OC * RT).
		oc 		: (float, int), optional
				Mask occultation.
				If present, value in range [0,1]
				If absent, 0.
		even_px
				: bool, optional
				If True, arrays are even-sized and the center of the
				telescope is taken at the 4 central pixels of the mask.
				Otherwise, the arrays are odd-sized and the center
				is at a certain pixel.
				This is ignored if mask width is defined, or the center
				is defined with (center_x, center_y).
		width 	: float, optional
				If width is > 0, the size of Rho, Phi, and Mask is defined.
				Otherwise, or if absent, the width is calculated automatically.
		center_x, center_y
				: float, optional
				If values are > 0, position of Rho, Phi, and Mask is defined.
				Otherwise, or if absent, then the center is calculated.
		verbose: bool (default = True)


	OUTPUTS:
		rho 	: (2D float numpy array) Radial coordinates
		phi 	: (2D float numpy array) Angular coordinates
		mask 	: (2D integer numpy array), 1 on the telescope and 0 elsewhere

	Note: 'useful_points' not implemented.
	"""




"""
Main Detection Algorithm function
"""
def detection_andromeda (snr_map, snr_norm = None, snr_stddev = None,
			flux_map = None, flux_stddev = None, flux_norm_stddev = None,
			threshold = None, size_subwindow = None,
			pixscale = None, oversampling = None,
			max_lobe_dist = None, owa = None, iwa = None,
			psf = None, img_int_time = None, psf_int_time = None,
			tnd = None, detection_files = None, nb_detections = None,
			verbose = True, obs_flag = 0, save_plots = False):
	

	if(verbose):
		print("\n\n ------ Andromeda Detection Algorithm ------ \n")
		print(" Input Checking...")

	"""
	------------------------------------------------------------------------------------------
	#00: INPUT CHECKING:
	------------------------------------------------------------------------------------------
	"""

	# Initialising maps class. Contains SNR and Flux maps.
	maps = Maps(snr_map, snr_norm, snr_stddev, flux_map, flux_stddev, flux_norm_stddev);
	maps.check_input(verbose)

	# Initialising data class. Contains all other important variables.
	data = Data(threshold, pixscale, oversampling, owa, iwa, size_subwindow, max_lobe_dist, obs_flag, maps.snr_map.shape[0])
	data.check_input(verbose, maps.flux_norm_stddev)

	if verbose:
		print("\n --- Info on parameters used ---")
		print(' Companions are sought from %.3f to %.3f lambda/D.' % (data.iwa, data.owa))
		print(' Threshold used for detection: %.3f.\n' % data.threshold)


	# Defining all variables stored in classes
	# Will continue cleaning up code afterwards
	sz = int(data.sz)
	size_snr = snr_map.shape[0]
	res = data.res
	iwa = data.iwa
	owa = data.owa
	max_lobe_dist = data.max_lobe_dist

	snr_stddev = np.ones_like(snr_map)
	snr_plot = np.copy(snr_norm)


	"""
	self.threshold = threshold
	self.pixscale = pixscale
	self.oversampling = oversampling
	self.owa = owa
	self.iwa = iwa
	# Size of subwindow from which to cut-out signals
	self.sz = size_subwindow
	# Size of input science images in pixels
	self.img_size = size_snr
	# Min distance between two signals to be considered separate
	self.max_lobe_dist = max_lobe_dist
	# Observation flag. Fit constraints for NaCo (0), NIRC2 (1), or STIM (2)
	self.obs_flag = obs_flag
	# Size of element of resolution in lambda/D
	self.res = 1.0
	# size of the subwindow in lambda/D
	self.sz_res = 1.0
	"""



	"""
	------------------------------------------------------------------------------------------
	#0: COMPANIONS' FIT PARAMETERS:
	------------------------------------------------------------------------------------------
	Constraints are set according to the sub-window size because 
	it is supposed to be chosen carefully by the user to fully enclose
	the planetary signal, as a function of the chosen delta_min parameter.
	"""

	#-Constraints for 2D-Gaussian fit:

	# Initialise class with "good" range of parameters
	# for a 2d gaussian fit
	gauss_params = Gauss_params(data.res, data.sz, data.obs_flag)

	#-Weight map for SNR 2D-Gaussian fit:
	sz_arr = idl.eclat( idl.distc(sz) )
	weight_snr = np.empty((sz,sz))

	weight_snr = 1 - sz_arr * (2./(sz-1))
	weight_snr = np.divide(weight_snr, np.amax(weight_snr))



	"""
	---------------------------
 	#1:	COMPANION DETECTION
	---------------------------
	"""

	# Redefining maps to avoid overwriting input images

	img = np.copy(maps.snr_map)
	norm = np.copy(maps.snr_norm)
	flux = np.copy(maps.flux_map)
	stddev_flux_norm = np.copy(maps.flux_norm_stddev)


	# Setting all negative image elements to zero for better fitting

	img[img < 0] = 0
	norm[norm < 0] = 0
	flux[flux < 0] = 0


	# Declaring arrays for detections data

	#Max detections in snr_map
	Pmax = 150

	snr_value = []
	nb_pixels = []
	flux_est = []
	flux_err = []
	flag_err_flux = []
	flag_err_pos = []
	# Error flags:
	# 0 - Good fit!
	# 2 - Fit failed
	# 3 - Too close to edge
	# 4 - Estimated local snr less than thresold

	# Detection oordinates in full map
	x_0 = []
	y_0 = []

  	#Detection coordinates in subwindow map, and their errors
	ax_0 = []
	ay_0 = []

	errx_0 = []
	erry_0 = []


	if np.amax(norm) < threshold:
		nb_detections = 0.
		print("Threshold is higher than max value of snr norm map. Decrease threshold!!")
		return

	# ---- All Detections Loop ----
	#Loops over every signal greater than threshold
	for r in range(0, Pmax):

		#Warnings
		if r == Pmax:
			print("WARNING: Too many companions found (", r+1, "). Increase Threshold!!")
			return

		if np.amax(norm) >= threshold:

			# Initialise position error flag
			flag_err_pos.append(0)

			# Calculating position indices of signal
			print("\n ---Signal above threshold found:", r+1)
			index_snr = idl.where(norm == np.amax(norm))
			index_x = index_snr % size_snr
			index_y = np.divide(index_snr, size_snr)
			i_ind = index_x[0]
			j_ind = index_y[0]

			# Calculating Mask
			rt = math.ceil(res)
			polaire_outputs = idl.polaire2(rt, width=size_snr, center_x=i_ind, center_y=j_ind)
			mask = polaire_outputs[2]

			# Get windows by cropping images with mask
			a1 = int(i_ind-(sz-1)/2.)
			a2 = int(i_ind+(sz-1)/2.)
			b1 = int(j_ind-(sz-1)/2.)
			b2 = int(j_ind+(sz-1)/2.)

			window_snr = (img * mask)[b1:b2,a1:a2]
			window_norm = 	(norm * mask)[b1:b2,a1:a2]
			window_flux = (flux * mask)[b1:b2,a1:a2] 
			window_stdev = (stddev_flux_norm * mask)[b1:b2,a1:a2]


			#----- Position Estimation -----

			# Test on the normal snr map
			ratio = snr_stddev[int(i_ind),int(j_ind)]
			ind = idl.where(window_snr > 0)
			ind2 = idl.where(window_snr >= threshold*ratio)

			# Number of pixels above threshold
			nb_pixels.append(ind2.shape[0])
			
			weight_here = weight_snr

			if (len(ind2.shape) == 1):
				weight_here = weight_here.flatten()
				for i in range (len(ind2)):
					weight_here[ind2[i]] = 1.
				reshape = weight_snr.shape[0]
				weight_here = weight_here.reshape((reshape,reshape))

			#--------- 2D Gaussian Fit

			gauss_snr = Gauss(window_snr)
			gauss_params.check(gauss_snr, verbose)

			gauss_norm = Gauss(window_norm)
			gauss_params.check(gauss_norm, verbose)


			#Plotting 3D fits of each detection
			if save_plots:
				plot_3d_window_gaussian(sz, gauss_snr, window_snr, filename=r)


			# outputs: gauss = (centroid_y, centroid_x, fwhm_y, fwhm_x, peak, theta)

			if ( int(gauss_snr.cx) >= window_snr.shape[0]) or ( int(gauss_snr.cy) >= window_snr.shape[1]):
				gauss_snr.good_fit = False

			if ( int(gauss_norm.cx) >= window_norm.shape[0]) or ( int(gauss_norm.cy) >= window_norm.shape[1]):
				gauss_norm.good_fit = False

			# ------ Switching between good/bad fits:

			errx_0.append(0.0)
			erry_0.append(0.0)

			# If SNR fit works:
			if (gauss_snr.good_fit):

				if verbose:
					print(" Planet candidate is fitted on the SNR Map")

				ax_0.append( int(gauss_snr.cx) )
				ay_0.append( int(gauss_snr.cy) )
				
				#
				# IDL mpfit2dpeak has more outputs than VIP equivalent!
				# Some code missing...
				#

				x_0.append( i_ind - (sz-1)/2. + gauss_snr.cx )
				y_0.append( j_ind - (sz-1)/2. + gauss_snr.cy )

				snr_value.append(window_norm[ax_0[r]][ay_0[r]])

				if snr_value[r] < threshold:
					flag_err_pos[r] = 4
				
				if verbose:
					print("Planet Candidate #", r+1)

				# Printing window_snr
				#plt.imshow(window_snr)
				#plt.show()

				dxmax = gauss_snr.fwhm_x * math.cos(-gauss_snr.angle)
				dymax = gauss_snr.fwhm_x * math.sin(-gauss_snr.angle)

				dxmin = -gauss_snr.fwhm_y * math.sin(-gauss_snr.angle)
				dymin = gauss_snr.fwhm_y * math.cos(-gauss_snr.angle)

			# If SNR Norm fit works
			elif (gauss_norm.good_fit):

				if verbose:
					print(" Planet candidate is fitted on the SNR Norm Map")

				ax_0.append( int(gauss_norm.cx) )
				ay_0.append( int(gauss_norm.cy) )
				
				#
				# IDL mpfit2dpeak has more outputs than VIP equivalent!
				# Some code missing...
				#

				x_0.append( i_ind - (sz-1)/2. + gauss_norm.cx )
				y_0.append( j_ind - (sz-1)/2. + gauss_norm.cy )

				ti = ax_0[r]

				tj = ay_0[r]

				print("wnorm = ", window_norm.shape)
				print(" ti=",ti," tj=",tj)

				tsnr_val = window_norm[ti][tj]

				snr_value.append( tsnr_val )

				print("Candidate position:", x_0[r], y_0[r])

				if snr_value[r] < threshold:
					flag_err_pos[r] = 4
				
				if verbose:
					print("Planet Candidate #", r+1)

				# Printing window_snr
				#plt.imshow(window_norm)
				#plt.show()

				dxmax = gauss_snr.fwhm_x * math.cos(-gauss_snr.angle)
				dymax = gauss_snr.fwhm_x * math.sin(-gauss_snr.angle)

				dxmin = -gauss_snr.fwhm_y * math.sin(-gauss_snr.angle)
				dymin = gauss_snr.fwhm_y * math.cos(-gauss_snr.angle)

			# If neither SNR nor Norm SNR map fits worked
			else:

				ax_0.append( int((sz-1)/2) )
				ay_0.append( int((sz-1)/2) )

				x_0.append( i_ind )
				y_0.append( j_ind )

				errx_0[r] = 0.5
				erry_0[r] = 0.5

				snr_value.append( window_snr[ax_0[r]][ay_0[r]] )

				flag_err_pos[r] = 2

			# Erasing companion to find the next one
			polaire_outputs = idl.polaire2(math.ceil(res), width=size_snr, center_x=x_0[r], center_y=y_0[r])
			mask_erase = polaire_outputs[2]
			
			norm *= (1 - mask_erase)
			img *= (1 - mask_erase)


			# Checking for companions too close to edge:
			dist = np.sqrt( np.abs( (size_snr-1)/2 - i_ind)**2 + abs( (size_snr-1)/2 - j_ind)**2)

			if dist > (2*owa*oversampling - 1*sz/2) or dist < (2*iwa*oversampling - 1*sz/3):
				if verbose:
					print("This companion is too close to the edge!")
				flag_err_pos[r] = 3

			# ---- Flux Estimation -----

			flux_est.append(0.0)
			flux_err.append(0.0)
			flag_err_flux.append(0)

			if flag_err_pos[r] == 0:

				gauss_flux = Gauss(window_flux)
				gauss_params.check(gauss_flux, verbose)

				if ( abs(sz-1)/2. - gauss_flux.cx <= gauss_params.cmax) or ( abs((sz-1)/2.- gauss_flux.cy) <= gauss_params.cmax ):

					print(" 2D Gaussian fit on flux pattern worked!")

					xp = (ax_0[r]-gauss_flux.cx)*math.cos(gauss_flux.angle) - (ay_0[r]-gauss_flux.cy)*math.sin(gauss_flux.angle)
					yp = (ax_0[r]-gauss_flux.cx)*math.sin(gauss_flux.angle) + (ay_0[r]-gauss_flux.cy)*math.cos(gauss_flux.angle)

					U = (xp / gauss_flux.fwhm_x)**2 + (yp/gauss_flux.fwhm_y)**2
					flux_est[r] = gauss_flux.peak * math.exp(-U/2.)
					# IDL: = constant + peak * exp(-U/2.)
					# No 'constant' return value from VIP gauss fit

					flag_err_flux[r] = 0
					flux_err[r] = 3.0*window_stdev[ idl.round(ax_0[r]) ][ idl.round(ay_0[r]) ]

				else:
					print(" 2D Gaussian fit on flux pattern did not work.")
					print(" Using rounded value instead.")
					flux_est[r] = window_flux[ idl.round(ax_0[r]) ][ idl.round(ay_0[r]) ]
					flag_err_flux[r] = 1
					flux_err[r] = 3.0*window_stdev[ idl.round(ax_0[r]) ][ idl.round(ay_0[r]) ]

			else:
				print(" No 2D Gaussian fit is performed in the flux map.")
				print(" Using rounded value instead.")
				flux_est[r] = window_flux[ idl.round(ax_0[r]) ][ idl.round(ay_0[r]) ]
				flag_err_flux[r] = 2
				flux_err[r] = 3.0*window_stdev[ idl.round(ax_0[r]) ][ idl.round(ay_0[r]) ]

		else:
			print("\n TOTAL NUMBER OF COMPANIONS FOUND: ", r, "\n")
			break
		# End of Main Loop


	#-- Recording results into arrays
	
	# Position
	Nb_planet = r
	coord_x = x_0
	coord_y = y_0
	errX = np.array(errx_0)
	errY = np.array(erry_0)

	# SNR
	snr_values = flux_est

	# Checking if companions are too close to each other
	# in which case, flag_pos = -1
	flag_pos = np.copy(flag_err_pos)

	for k in range(0, Nb_planet):
		for h in range(0, Nb_planet):
			if k==h:
				continue
			x_dist = coord_x[k] - coord_x[h]
			y_dist = coord_y[k] - coord_y[h]
			distance = math.sqrt( (x_dist)**2 + (y_dist)**2 )

			if (distance < max_lobe_dist):
				# Signals too close to each other!!
				# Set falg of weaker one to (-1)
				if (snr_value[k] > snr_value[h]):
					flag_pos[h] = -1
				else:
					flag_pos[k] = -1


	#-- PLOTTING DETECTION MAP

	plt.imshow(snr_plot)
	plt.title('Detections in SNR Map')
	
	g_patch = mpatches.Patch( color='w', label='Good Fit')
	k_patch = mpatches.Patch( color='c', label='Too close to other signal')
	y_patch = mpatches.Patch( color='y', label='Too close to edge')
	r_patch = mpatches.Patch( color='r', label='Bad Fit')
	w_patch = mpatches.Patch( color='k', label='Low SNR')

	plt.legend(handles=[g_patch, k_patch, y_patch, r_patch, w_patch], 
				loc = 'upper center',
				bbox_to_anchor=(1.05, 1),
				prop = {'size': 8},
				ncol = 1, fancybox = True, shadow = True)

	for i in range(0, len(ax_0)):
		color = 'w'
		# No errors
		if flag_pos[i] == 0:
			color = 'w'
		# Bad fit
		elif flag_pos[i] == 2:
			color = 'r'
		# Too close to edge
		elif flag_pos[i] == 3:
			color = 'y'
		# Good fit
		elif flag_pos[i] == 4:
			color = 'k'
		# Too close to other signal
		elif flag_pos[i] == -1:
			color = 'c'

		circ = plt.Circle((coord_x[i],coord_y[i]),10, fill=False, color=color)
		plt.gcf().gca().add_artist(circ)

	plt.show()
	#fig.savefig("detections.png")


	#--------------------------------
	# 3 - DISPLAY RESULTS: POSITION
	#--------------------------------

	# 1 - Entire List (using Pandas Dataframe table)
	print("\n")
	print(" POSITION INFORMATION:")

	d = {'x': coord_x, 'y': coord_y, 'snr': snr_value, 'pixels-nb': nb_pixels, 'flag': flag_pos}

	df = pd.DataFrame(data=d)
	print(df)

	df.to_csv("detections.csv")


	#--------------------------------
	# 4 - DISPLAY RESULTS: FLUX
	#--------------------------------

	# Position
	offset_x = np.empty(Nb_planet, dtype='float')
	offset_y = np.empty(Nb_planet, dtype='float')

	# Angular Separation
	pos_arcsec = np.empty(Nb_planet, dtype='float')
	errpos_arcsec_pos = np.empty(Nb_planet, dtype='float')
	errpos_arcsec_neg = np.empty(Nb_planet, dtype='float')
	pos_arcsec_err = np.empty(Nb_planet, dtype='float')

	# Position angle
	pos_angle = np.empty(Nb_planet, dtype='float')
	errpos_angle_pos = np.empty(Nb_planet, dtype='float')
	errposs_angle_neg = np.empty(Nb_planet, dtype='float')
	pos_angle_err = np.empty(Nb_planet, dtype='float')

	# Contrast and magnitude
	contrast = np.empty(Nb_planet, dtype='float')
	contrast_err = np.empty(Nb_planet, dtype='float')
	delta_lam = np.empty(Nb_planet, dtype='float')
	delta_lam_err = np.empty(Nb_planet, dtype='float')


	# ----------Calculating

	coord_x = np.array(coord_x)
	coord_y = np.array(coord_y)

	# Angular Separations
	align = (size_snr-1)/2.0
	offset_x = align - coord_x
	offset_y = align - coord_y

	pos_arcsec = pixscale * np.sqrt( (-offset_x)**2 + (-offset_y)**2 )

	# Position Angles
	factors = np.ones(Nb_planet)
	factors[ coord_x > align ] = 3
	pos_angle = ( np.arctan2(offset_y, offset_x) + factors*np.pi/2.) * 180.0/np.pi

	# Position Errors
	sepx = np.empty(Nb_planet)
	sepy = np.empty(Nb_planet)

	sepx[ coord_x > align ] = -offset_x[ coord_x > align ]
	sepx[ coord_x <= align] = offset_x[ coord_x <= align]

	sepy[ coord_y > align ] = -offset_y[ coord_y > align ]
	sepy[ coord_y <= align] = offset_y[ coord_y <= align]

	sep_all = np.sqrt( sepx**2 + sepy**2 )
	pos_arsec_err = pixscale * (1/sep_all) * (sepx*errX + sepy*errY)

	err_pa = (sepx*sepy) / (sepx**2+sepy**2) * (-errX/sepx + errY/sepy)
	pos_angle_err = np.fabs( err_pa * 180.0 / np.pi )

	contrast = np.copy(flux_est)
	contrast_err = np.copy(flux_err)

	# Contrast in terms of magnitude
	delta_lam = np.fabs( 2.5*np.log10(contrast) )

	# Contrast Errors
	contrast_max = contrast + contrast_err
	contrast_min = contrast - contrast_err

	delta_lam_max = np.fabs(2.5*np.log10(contrast_max))
	delta_lam_min = np.fabs(2.5*np.log10(contrast_min))

	for p in range(Nb_planet):
		delta_lam_err[p] = max(delta_lam_min[p] - delta_lam[p], delta_lam[p] - delta_lam_max[p])

	#min_condition = delta_lam_min - delta_lam
	#max_condition = delta_lam - delta_lam_max

	#delta_lam_err[ min_condition > max_condition ] = min_condition
	#delta_lam_err[ max_condition >= min_condition ] = max_condition


	# ---------- Printing Results
	d = {'arcsec': pos_arcsec, 'angle': pos_angle, 'contrast': contrast, 'mag contrast': delta_lam, 'flag': flag_err_flux}
	df = pd.DataFrame(data=d)
	print(df)
	df.to_csv("fluxes.csv")

	print("DONE")
	return




def flip(data):
	data = np.flip(data, axis=0)
	return data



def main():

	print("\n === Andromeda Detection Algorithm ===")

	# Ask for image type to analyse
	option = 0
	while(True):
		print("\n Analyse...")
		print(" 1) Andromeda outputs")
		print(" 2) A specific image")
		option_str = input("")
		try:
			option = int(option_str)
		except ValueError:
			print(" Error: invalid option")
		if (option == 1 or option == 2):
			break

	# Ask for images directory
	path = ''
	starname = ''
	starpath = ''

	snrmap_path = ''
	snrnorm_path = ''

	while(not os.path.isfile(starpath) and option == 2):
		starpath = input(" Input image directory: ")
		if (not os.path.isfile(starpath)):
			print(" File '%s' not found" % starpath)

	while(not os.path.isfile(starpath) and option == 1):
		path = input(" Input folder name where data is located (example: 'hr_8799' or '.\\'): ")
		if(path[-1] != '\\') or (path[-1] != '/'):
			path += '\\'
		starname = input(" Input star name: ")
		
		#Check that SNR maps exist
		snrmap_path = path + "andr_snrmap_" + starname + ".fits"
		snrnorm_path = path + "andr_snrnorm_" + starname + ".fits"

		if(not os.path.isdir(path)):
			print(" Error: path '%s' does not exist" % path)
			continue

		if (not os.path.isfile(snrmap_path) and not os.path.isfile(snrnorm_path)):
			print(" Error: SNR and SNR-norm maps not found:")
			print(" ->", snrmap_path)
			print(" ->", snrnorm_path)
			continue

		break

	snr_map = None
	snr_norm = None
	flux_map = None
	stddev_flux_map = None
	stddev_flux_norm = None

	if (option == 1):
		flux_map_path = path + "andr_contrast_" + starname + ".fits"
		stddev_flux_map_path = path + "andr_contrast_stddev_" + starname + ".fits"
		stddev_flux_norm_path = path + "andr_contrastnorm_stddev_" + starname + ".fits"

		if (os.path.isfile(snrmap_path)):
			snr_map = np.flip(vip.fits.fits.open_fits(snrmap_path), axis=0)

		if (os.path.isfile(snrnorm_path)):
			snr_norm = np.flip(vip.fits.fits.open_fits(snrnorm_path), axis=0)

		if (os.path.isfile(flux_map_path)):
			flux_map = np.flip(vip.fits.fits.open_fits(flux_map_path), axis=0)

		if (os.path.isfile(stddev_flux_map_path)):
			stddev_flux_map = np.flip(vip.fits.fits.open_fits(stddev_flux_map_path), axis=0)

		if (os.path.isfile(stddev_flux_norm_path)):
			stddev_flux_norm = np.flip(vip.fits.fits.open_fits(stddev_flux_norm_path), axis=0)

		if (not os.path.isfile(snrmap_path)):
			snr_map = np.copy(snr_norm)

	elif (option == 2):
		snr_map = np.flip(vip.fits.fits.open_fits(starpath), axis=0)
		snr_norm = np.copy(snr_map)
		flux_map = np.ones(snr_map.shape)
		stddev_flux_map = np.ones(snr_map.shape)
		stddev_flux_norm = np.ones(snr_map.shape)

	# Pixscale for K-band
	Kp_wavelength = 2.124E-6
	pixscale_keck = 0.00953
	diam_tel = 10
	pixscale = pixscale_keck * 1E3
	pixscale_nyq = 1/2 * (180*3600*1E3 / np.pi) * (Kp_wavelength / diam_tel)
	oversampling = pixscale_nyq / pixscale
	print(" Oversampling = ",oversampling)

	# ---Threshold
	# For Naco/VLT ~= 5
	# For Nirc2/Keck ~= 10
	threshold = float(0.0)
	while(True):
		th = input(" Input threshold: ")
		try:
			th = float(th)
		except ValueError:
			print(" Error: threshold must be int or float")
			continue

		if (th <= 0):
			print(" Error: threshold must be greater than zero")
			continue
		threshold = th
		break

	# ---Flag for constraints
	# 0 - NaCo-VLT
	# 1 - NIRC2-Keck
	# 2 - STIM on NIRC2-Keck
	obs_flag = 1

	detection_andromeda (snr_map, snr_norm = snr_norm, snr_stddev = None,
			flux_map = flux_map, flux_stddev = stddev_flux_map, flux_norm_stddev = stddev_flux_norm,
			threshold = threshold, size_subwindow = None,
			pixscale = pixscale, oversampling = oversampling,
			max_lobe_dist = None, owa = None, iwa = None,
			psf = None, img_int_time = None, psf_int_time = None,
			tnd = None, verbose = True, obs_flag = obs_flag, save_plots = False)



if __name__ == '__main__':
	main()
	exit()
