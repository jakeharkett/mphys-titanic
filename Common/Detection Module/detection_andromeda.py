import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import math
import numpy as np
import pandas as pd
import os
import sys
import argparse

import vip_hci as vip

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib import cm


"""

	ANDROMEDA DETECTION ALGORITHM

	Algorithm for automatic detection of planetary signals
	in high-contrast imaging datasets.
	
	Developed by Faustine Cantalloube in IDL.
	Translated to Python by Jorge Fernandez (jorgefz.fernandez@gmail.com)

	Implemented from version:
		Revision: 1.5, Date: 2018/06/25 

	Check the paper Cantalloube et al. (2015) for more info:
	https://www.aanda.org/articles/aa/abs/2015/10/aa25571-14/aa25571-14.html


	RUNNING THE PROGRAM

	For default options, simply run the script in Python
	with one of the following commands:
		py detection_andromeda.py
		python3 detection_andromeda.py

	Then, you will be prompted to schoose between analysing
	an Andromeda Output set (1) or a single image (2).

	Choosing to analyse the Andromeda Output, you will be asked
	to input the folder in which the images can be found.
	If the images are in the current forlder, simply type a dot '.'
	Then, input the starname as seen in the files' names.
	e.g.	ANDR_snrmap_HR_8799.fits has star name 'HR_8799'.

	Choosing to analyse a single image will simply prompt you
	to type the path to the image.

	Finally, you will be asked for a threshold. Only signals
	with a Signal-to-Noise ratio larger than this number will be
	sought in the Normalized SNR map (or the regular SNR map
	if the latter is missing).

	Running the script will generate a map 'detections.png' wherein
	the signals are circled and colour-coded, 
	as well as a file 'detections.csv' with signal data 
	and 'fluxes.csv' with their flux data
	(if flux map was provided).
	These files will be generated on the path provided, where
	the input dataset or image is located.


	COMMAND LINE ARGUMENTS

	The script may be called with a number of command line arguments.
	To check the available commands, run either of the following:
		py detection_andromeda.py -h
		py detection_andromeda.py --help

	The commands are the following:

		usage: newdet.py [-h] [-v] [-p] [-g GAUSS] [-m MAXDET]

		optional arguments:
		  -h, --help            show this help message and exit
		  -v, --verbose         Displays extra information during companion detection.
		                        Default = False
		  -p, --plot3d          Generates 3D plots of signal and gaussian fits.
		                        Default = False
		  -g GAUSS, --gauss GAUSS
		                        Specifies set of gaussian constraints to use on
		                        signals.Names: naco, nirc2, stim, other
		  -m MAXDET, --maxdet MAXDET
		                        Max limit of detections in the map. Default = 100

	Examples:

		To run the script with verbose messages (extra info):
			py detection_andromeda.py -v
			py detection_andromeda.py --verbose

		To run the script with the verbose flag and save 3D plots:
			py detection_andromeda.py -v -p
			py detection_andromeda.py --verbose --plot3d

		To run the script with a constraint set 'naco':
			py detection_andromeda.py -g naco
			py detection_andromeda.py --gauss naco

		To run the script with 50 max detections:
			py detection_andromeda.py -m 50
			py detection_andromeda.py --maxdet 50


	FURTHER DEVELOPMENT

		- Improve Gaussian Constraints for Keck data
		- Improve Gaussian Constraints for input STIM maps

"""

# Super useful video -> https://www.youtube.com/watch?v=dQw4w9WgXcQ

"""
Stores the data maps used in detection
"""
class MapsClass:
	# Input data maps
	def __init__(self, smap=None, snorm=None, sdev=None,
					fmap=None, fdev=None, fndev=None, verbose=False):
		# SNR maps
		self.smap = smap		# SNR map
		self.snorm = snorm		# SNR normalized map
		self.sdev = sdev		# SNR standard deviation map
		# Flux maps
		self.fmap = fmap		# Flux map		
		self.fdev = fdev 		# Flux standard deviation map
		self.fndev = fndev		# FLux normalized standard deviation map

		# Flags
		self.noflux = False		#Indicates if flux map is present

		# Map checks
		if not isNumpyArray(self.smap) and not isNumpyArray(self.snorm):
			print(" Error: both SNR and Normalized SNR map are missing")
			exit(-1)
			#raise ValueError(" Error: both SNR and Normalized SNR map are missing")
		elif isNumpyArray(self.smap) and not isNumpyArray(self.snorm):
			self.snorm = np.copy(self.smap)
		elif isNumpyArray(self.snorm) and not isNumpyArray(self.snorm):
			self.smap = np.copy(self.snorm)

		size = self.smap.shape	#Temp variable, size of images

		if not isNumpyArray(self.sdev):
			self.sdev = np.zeros( size, dtype=float )
			if(verbose):
				print(" Missing SNR standard deviation map. Ignoring...")

		if not isNumpyArray(self.fmap):
			self.noflux = True
			if(verbose):
				print(" Missing flux map. No flux calculations will be performed.")
		else:
			if not isNumpyArray(self.fdev):
				self.fdev = np.zeros( size, dtype=float )
				if(verbose):
					print(" Missing standard deviation flux map. Ignoring...")
			if not isNumpyArray(self.fndev):
				self.fndev = np.zeros( size, dtype=float )
				if(verbose):
					print(" Missing normalized standard deviation flux map. Ignoring...")

		# Make all negative values in maps zero
		self.smap[self.smap < 0] = 0
		self.snorm[self.snorm < 0] = 0
		if(not self.noflux):
			self.fmap[self.fmap < 0] = 0


"""
Stores physical variables about data
"""
class VarsClass:
	def __init__(self, limit=0, pxscale=0, oversamp=0,
				owa=0, iwa=0, subw_sz=0, map_sz=0, neigh=0,
				verbose=False, plot3d=False, path="./"):

		self.limit = limit			# Threshold for detections in map
		self.pxscale = pxscale		# Pixel scale of data
		self.oversamp = oversamp	# Oversampling factor
		self.owa = owa				# Outer Working Angle
		self.iwa = iwa				# Inner Working Angle

		self.subw_sz = subw_sz		# Size of signal-cropping subwindow
		self.map_sz = map_sz		# Size of input images
		self.res = 0 				# Size of element of resolution in lambda/D
		self.subw_res = 0			# Sizeo of subwindow in lambda/D

		self.neigh = neigh		# Max distance to neighbour signals
								# Closer and the weaker signal is ignored

		self.plot3d = plot3d 	# Flag to generate 3D plots of Signal and Gaussian fit
		self.verbose = verbose	# Flag to show extra info during detection
		self.path = path		# Path to star data


		#Check input variables
		if( limit <= 0):
			print(" Error: threshold must be greater than zero")
			exit(-1)

		if( pxscale <= 0):
			print(" Error: pixel scale must be greater than zero")
			exit(-1)

		if( oversamp <= 0):
			print(" Error: oversampling factor must be greater than zero")
			exit(-1)

		if(self.verbose):
			print(" Oversampling set to",self.oversamp, "lambda/D")

		if( map_sz <= 0):
			print(" Error: map size must be greater than zero")
			exit(-1)

		if( owa <= 0 ):
			self.owa = (self.map_sz/2. - 1) / (2.*self.oversamp)
			if(self.verbose):
				print(" OWA set to", self.owa, "lambda/D")
		
		if( iwa <= 0):
			self.iwa = 1.0
			if(self.verbose):
				print(" IWA set to", self.iwa, "lambda/D")

		if( subw_sz <= 0):
			self.subw_sz = math.floor(4*self.oversamp) + 2
			if(self.verbose):
				print(" Subwindow size set to", self.subw_sz, "pixels")

		self.res = 2*self.oversamp
		self.subw_res = float(self.subw_sz)/float(self.res)

		if( neigh <= 0):
			self.neigh = math.floor(10*self.oversamp)
			if(self.verbose):
				print(" Max neighbour distance set to", self.neigh, "lambda/D")


"""
Stores 'good' Gaussian fit constraints
"""
class GaussParamsClass:
	def __init__(self, res, subw_sz, telescope_type):
		self.sz = subw_sz
		self.lsize = [0,0]	# [Max,Min] range for larger FWHM
		self.ssize = [0,0]	# [Max,Min] range for smaller FWHM
		self.cmax = 0		# Max deviation from subwindow centre
		self.fwhm_ratio = 0	# Max ratio between Large and Small FWHMs
		self.fwhm_min = 0	# Min value for small FWHM
		self.fwhm_max = 0	# Max value for large FWHM
		self.ttype = telescope_type

		"""
		FWHM = Full Width Half Maximum
		In 2D Gaussian functions, there are two FWHMs, the X and Y ones.
		"""

		""" 
		'telescope_type' defines constraints type
		e.g. contraints for Keck-NIRC2 are different 
		from those of VLT-NaCo.
		Feel free to edit these values or add new ones
		to improve Goodness-of-Gaussian-fit accuracy
		"""
		if( telescope_type == "nirc2" ):
			# Max,Min limits for larger FWHM
			self.lsize = [ 1.10*res, 0.25*res ]
			# Max,Min limits for smaller FWHM
			self.ssize = [ 0.90*res, 0.15*res ]
			# Max deviation of fit centre from subwindow centre
			self.cmax = 1.50
			# Max fwhm large / small ratio (max signal enlargement)
			self.fwhm_ratio = 2.
			# Min and max values for FWHMs, secondary check
			self.fwhm_min = (0.2)*subw_sz
			self.fwhm_max = 1.5*subw_sz

		elif( telescope_type == "naco" ):
			self.lsize = [ 0.95*res, 0.25*res ]
			self.ssize = [ 0.65*res, 0.20*res ]
			self.cmax = 0.75
			#-----------------------
			self.fwhm_ratio = 2.
			self.fwhm_min = (0.2)*subw_sz	
			self.fwhm_max = 1.5*subw_sz

		elif( telescope_type == "stim" ):
			self.lsize = [ 0.60*res, 0.25*res ]
			self.ssize = [ 0.65*res, 0.15*res ]
			self.cmax = 1.50
			self.fwhm_ratio = 2.
			self.fwhm_min = (0.2)*subw_sz
			self.fwhm_max = 1.5*subw_sz

		elif( telescope_type == "other"):
			self.lsize = [ 0.60*res, 0.25*res ]
			self.ssize = [ 0.65*res, 0.15*res ]
			self.cmax = 1.50
			self.fwhm_ratio = 2.
			self.fwhm_min = (0.2)*subw_sz
			self.fwhm_max = 1.5*subw_sz

		else:
			print(" Error: unknown telescope type: '%s'" % telescope_type)
			exit(-1)

	def check_fit(self, fwhm, center, subw_sz):
		# Find larger/smaller between FWHM_x and FWHM_y
		fwhm_large = 0
		fwhm_small = 0
		if(fwhm[0] >= fwhm[1]):
			fwhm_large = fwhm[0]
			fwhm_small = fwhm[1]
		else:
			fwhm_large = fwhm[1]
			fwhm_small = fwhm[0]

		# Compare with constraints
		subw_center = (subw_sz-1)/2.0

		if (fwhm_small >= self.ssize[0]) or (fwhm_small <= self.ssize[1]):
			return False

		if (fwhm_large >= self.lsize[0]) or (fwhm_large <= self.lsize[1]):
			return False
		
		if ( abs(subw_center-center[0]) > self.cmax ):
			return False

		if ( abs(subw_center-center[1]) > self.cmax ):
			return False

		return True

	def check_SNR_fit(self, Signal):
		subw_sz = Signal.maps.smap.shape[0]
		fwhm = [Signal.fwhm_x, Signal.fwhm_y]
		center = [Signal.cx, Signal.cy]
		Signal.goodfit = self.check_fit(fwhm, center, subw_sz)

	def check_flux_fit(self, Signal):
		subw_sz = Signal.maps.fmap.shape[0]
		fwhm = [Signal.f_fwhm_x, Signal.f_fwhm_y]
		center = [Signal.f_cx, Signal.f_cy]
		Signal.f_goodfit = self.check_fit(fwhm, center, subw_sz)
		


class DetDataClass:
	def __init__(self, map_class, var_class, gausspar_class):
		self.Map = map_class
		self.Var = var_class
		self.GaussParam = gausspar_class
		self.SignalArray = []

		
"""
Stores data of a detected signal
"""
class SignalClass:
	def __init__(self, wmap=None, wnorm=None, wdev=None,
					wflux=None, wfdev=None, wfndev=None):
		# Maps class stores subwindow on different maps
		self.maps = MapsClass(smap=wmap, snorm=wnorm, sdev=wdev,
					fmap=wflux, fdev=wfdev, fndev=wfndev,
					verbose=False)

		self.pos = [0,0]	# Position in full map
		self.ind_max = [0,0] 	# Indices of maximum in full map

		# SNR MAP - Gaussian fit parameters for signal
		self.snr_val = 0	# Signal SNR value
		self.snr_err = 0

		self.goodfit = False
		self.cx = 0
		self.cy = 0
		self.fwhm_x = 0
		self.fwhm_y = 0
		self.peak = 0
		self.angle = 0

		self.pos_flag = 0
			# -1 - Good/Bad fit on SNR/flux
			# 0 - Good fit on both
			# 1 - Bad fit on both
			# 2 - Close to edge
			# 3 - Close to neighbour
			# 4 - Low SNR

		# FLUX MAP - Gaussian fit parameters for signal
		self.flux_val = 0	# Signal flux value
		self.flux_err = 0

		self.f_goodfit = False
		self.f_cx = 0
		self.f_cy = 0
		self.f_fwhm_x = 0
		self.f_fwhm_y = 0
		self.f_peak = 0
		self.f_angle = 0

		# CONTRAST
		self.pos_arcsec = 0
		self.pos_arcsec_err = 0

		self.pos_angle = 0
		self.pos_angle_err = 0

		self.contrast = 0
		self.contrast_err = 0


	# Fits 2D gaussian to the signal and tests constraints
	def fit_gaussian(self, DetData, detmap):
		gauss_out = vip.var.fit_2dgaussian( self.maps.smap,
					full_output = True, debug=False)
		# Centroid
		self.cx = gauss_out.iloc[0,1]
		self.cy = gauss_out.iloc[0,0]
		
		#Full Width Half Maxima
		self.fwhm_x = gauss_out.iloc[0,3]
		self.fwhm_y = gauss_out.iloc[0,2]

		self.peak = gauss_out.iloc[0,4]
		self.angle = gauss_out.iloc[0,5]

		#Check with constraints
		DetData.GaussParam.check_SNR_fit(self)
		if(not self.goodfit) and (DetData.Var.verbose):
			print(" 	> Bad Gaussian fit on SNR Map!")

		# Set peak position in full map
		subw_sz = self.maps.smap.shape[0]
		self.ind_max = get_max_ind(detmap)
		pos_x = self.ind_max[0] - (subw_sz-1)/2.0 + self.cx
		pos_y = self.ind_max[1] - (subw_sz-1)/2.0 + self.cy
		self.pos = [int(pos_x), int(pos_y)]

		# Set SNR value
		self.snr_val = DetData.Map.snorm[self.pos[1]][self.pos[0]]
		self.snr_err = DetData.Map.sdev[self.pos[1]][self.pos[0]]
		if DetData.Var.verbose:
			print(" 	> SNR of signal is", self.snr_val)

	# Sets signal positional flag
	# Only to be called when detection ends
	def set_pos_flag(self, DetData, ind):
		# Good SNR / Bad FLux fits, or vice-versa
		if(self.goodfit and not self.f_goodfit) or \
			(not self.goodfit and self.f_goodfit):
			self.pos_flag = -1
		# Good Fit on both
		elif(self.goodfit and self.f_goodfit):
			self.pos_flag = 0
		# Close to edge
		if( is_close_to_edge(DetData, ind) ):
			self.pos_flag = 2
		# Close to neigh
		elif( is_close_to_neighbour(DetData, ind) ):
			self.pos_flag = 3
		# Low SNR
		elif( self.snr_val < DetData.Var.limit ):
			self.pos_flag = 4
		# Bad Fit on both
		elif(not self.goodfit and not self.f_goodfit):
			self.pos_flag = 1

	def fit_gaussian_flux(self, DetData):
		if( self.maps.noflux ):
			return

		gauss_out = vip.var.fit_2dgaussian( self.maps.fmap,
					full_output = True, debug=False)
		# Centroid
		self.f_cx = gauss_out.iloc[0,1]
		self.f_cy = gauss_out.iloc[0,0]
		
		#Full Width Half Maxima
		self.f_fwhm_x = gauss_out.iloc[0,3]
		self.f_fwhm_y = gauss_out.iloc[0,2]

		self.f_peak = gauss_out.iloc[0,4]
		self.f_angle = gauss_out.iloc[0,5]

		#Check with constraints
		DetData.GaussParam.check_flux_fit(self)
		if(not self.f_goodfit) and (DetData.Var.verbose):
			print(" 	> Bad Gaussian fit on Flux Map!")

		if( self.f_goodfit):
			# Calculate flux from flux value at Gaussian distribution center
			xdiff = self.cx - self.f_cx
			ydiff = self.cy - self.f_cy
			xp = xdiff*math.cos(self.f_angle) - ydiff*math.sin(self.f_angle)
			yp = xdiff*math.cos(self.f_angle) + ydiff*math.sin(self.f_angle)
			expn = (xp/self.f_fwhm_x)**2 + (yp/self.f_fwhm_y)**2
			self.flux_val = self.f_peak * math.exp(expn/2.)

		else:
			# Otherwise, estimate flux from subwindow
			self.flux_val = DetData.Map.fmap[ self.pos[1] ][ self.pos[0] ]
		
		self.flux_err = 3.0*DetData.Map.fndev[ self.pos[1] ][ self.pos[0] ]
		if DetData.Var.verbose:
			print(" 	> Flux of signal is", self.flux_val)


		# CONTRAST CALCULATIONS
		pxscale = DetData.Var.pxscale
		subw_cent = (DetData.Var.subw_sz - 1)/2.0
		offset_x = self.pos[0] - subw_cent
		offset_y = self.pos[1] - subw_cent

		self.pos_arcsec = pxscale * math.sqrt( offset_x**2 + offset_y**2 )

		factor = 1
		if (self.pos[0] > subw_cent):
			factor = 3

		self.pos_angle = (np.arctan2(offset_y, offset_x) + factor*np.pi/2.) * 180.0/np.pi

		sepx = abs(offset_x)
		sepy = abs(offset_y)
		sepall = math.sqrt(sepx**2 + sepy**2)
		# errX and errY come from Gaussian Fit
		errX = 0
		errY = 0

		self.pos_arcsec_err = pxscale * ( 1/sepall *(sepx*errX + sepy*errY) )
		
		errpa = (sepx*sepy)/(sepx**2+sepy**2) * (-errX/sepx + errY/sepy)
		self.pos_angle_err = abs(errpa * 180 / math.pi)

		# Flux errors
		flux_max = self.flux_val + self.flux_err
		flux_min = self.flux_val - self.flux_err

		if(flux_min <= 0):
			flux_min = self.flux_val
			if(DetData.Var.verbose):
				print(" Warning: flux error is too large. Ignoring...")

		self.contrast = np.fabs( 2.5*np.log10(self.flux_val))

		contrast_max = np.fabs( 2.5*np.log10(flux_max))
		contrast_min = np.fabs( 2.5*np.log10(flux_min))
		self.contrast_err = max(abs(self.contrast-contrast_min), abs(self.contrast - contrast_max))

		if DetData.Var.verbose:
			print(" 	> Contrast of signal is", self.contrast)





# =================================================================
#		FUNCTIONS
# =================================================================


# ============= IDL functions ================= 
class idl:

	"""
	Implementation of POLAIRE2 function from IDL. Author: L. Mugnier
	Generates a circular mask using a telescope radius.
	PARAMETERS
		rt 	: Positive float or int. Telescope radius in pixels.  
		leq	: Bool. If True, points at distance of telescope radius are included.
			-> True = Mask defined by points of (val <= RT) to (val >= OC*RT)
			-> False = Mask defined by points of (val < RT) to (val >= OC*RT)
		oc 	: Float or int in range [0,1]. Mask occultation percentage.
		width 	: Int. Size of mask in pixels. Can be calculated automatically.
		even_px	: Bool. Defines wether mask is even or odd.
			-> True = Centre of telescope is defined as 4 central pixels of mask.
			-> False = Centre of telescope is simply central pixel.
		center_x	: Float.
		center_y	: Float. Manually set centre x,y of mask
		verbose 	: Bool. If True, print information during calculation.

	OUTPUTS:
		rho 	: 2D float numpy array. Radial coordinates
		phi 	: 2D float numpy array. Angular coordinates
		mask 	: 2D integer numpy array. Has values 1 on the telescope aperture and 0 elsewhere
	"""
	def polaire2(rt = 0, width = 0, oc = 0.0,
				center_x = None, center_y = None,
				even_px = False, leq = False,
				verbose = True):

		# 	Input Checking
		if (rt <= 0):
			print(" POLAIRE: Telescope radius must be greater than zero")
			return

		if (oc < 0) or (oc > 1):
			print(" POLAIRE: Occultation must be in range [0,1]. Will be set to 0")
			oc = 0.0

		# Checking Width
		if (width <= 0):
			if even_px:
				width = int(2 * idl.round(rt))
			else:
				width = int(2 * int(rt) + 1)
			if(verbose):
				print(" POLAIRE: Width set to", width, "pixels")

		# Definiton of center x,y
		if not center_x:
			center_x = float( (width-1.0)/2 )
		if not center_y:
			center_y = float( (width-1.0)/2 )
		if verbose:
			print(" POLAIRE: Center_xy set to", center_x, center_y)


		# Calculating Angles Rho and Phi
		x = np.arange(width*width, dtype=float)
		x = x.reshape((width,width))
		x = x % width
		# These 3 lines are equivalent to IDL:
		#	x = float( dindgen(widht,width) mod width )

		y = np.transpose(x)
		x = x - center_x
		y = y - center_y

		rho = np.sqrt(x**2 + y**2)/float(rt)
		phi = np.arctan2(y, x + (rho==0).astype(int) )

		# Calculating Mask
		if leq:
			mask = ( (rho <= 1)*(rho >= oc) ).astype(int)
		else:
			mask = ( (rho < 1)*(rho >= oc) ).astype(int)

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
				if (x[i] >= (n-x)[i]):
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

# ========== End of IDL Functions =====================



# Checks if input variable is a Numpy Array. Returns Bool.
def isNumpyArray(arr):
	return isinstance(arr, np.ndarray)

# Creates a directory if it didn't exist before
def makeDir(path):
	if( not os.path.isdir(path) ):
		os.mkdir(path)


"""
Makes a 3D plot with a data image and its gaussian fit
PARAMETERS:
	sz 		: Size of signal
	fwhm 	: Array of lenth 2. X-Y values for FWHM of Gaussian fit
	signal 	: Signal data array
	filename: Output filename for image
"""
def plot_3d_window_gaussian(sz=0, fwhm=None, window=None, 
							path=None, filename=None):

	# Checking filename input
	if (filename == None):
		filename = "signal"

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
	ax.title.set_text('Signal')
	fig.colorbar(surf, shrink=0.5, aspect=10)

	# ----------------
	# 2nd Plot: Gaussian Fit
	# ----------------
	ax = fig.add_subplot(1, 2, 2, projection='3d')
	ax.title.set_text('Gaussian Fit')
	
	# Calculate gaussian parameter sigma
	sigma_x = fwhm[0] / ( 2*math.sqrt(2*math.log(2)) )
	sigma_y = fwhm[1] / ( 2*math.sqrt(2*math.log(2)) )
	# Shift center
	sh = (sz-1)/2.
	# Calculate gaussian curve values on meshgrid
	fit = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((x-sh)**2/(2*sigma_x**2)
    	+ (y-sh)**2/(2*sigma_y**2))))

	surf2 = ax.plot_surface(x, y, fit, cmap='Blues')
	fig.colorbar(surf2)

	makeDir(path + "detections")
	fig.savefig(path + "detections/"+filename+".png")
	#plt.show()
	plt.close(fig)


# Retrieves xy tuple with indices of maximum in map
def get_max_ind(_map):
	flat_index_max = np.argmax(_map, axis=None)
	index_max = np.unravel_index(flat_index_max, _map.shape)
	index_max = np.flip(index_max)
	return index_max

# Given
def crop_window(_map, mask, wx, wy):
	# Apply mask
	masked_map = np.copy(_map) * mask
	# Crop out subwindow using input indices
	window = masked_map[ wy[0]:wy[1], wx[0]:wx[1] ]
	return window

# Determines whether signal is too close to edge of full map
def is_close_to_edge(DetData, ind):
	currSignal = DetData.SignalArray[ind]
	size_snr = DetData.Map.smap.shape[0]
	center = (size_snr-1)/2.
	pos_x = center - currSignal.ind_max[0]
	pos_y = center - currSignal.ind_max[1]

	dist = np.sqrt( pos_x**2 + pos_y**2 )

	bmin = 2*DetData.Var.iwa*DetData.Var.oversamp - 1*DetData.Var.subw_sz/3
	bmax = 2*DetData.Var.owa*DetData.Var.oversamp - 1*DetData.Var.subw_sz/2.

	if (dist > bmax) or (dist < bmin):
		return True
	return False
	
# Calculates 2D distance between two points
def get_distance(a, b):
	dist = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)
	return dist

# Determines whether signal is too close to another one
# and which one is the stronger one.
def is_close_to_neighbour(DetData, ind):

	currSignal = DetData.SignalArray[ind]
	neighSignals = []

	# Loop over all signals and pick those within limit
	for i,s in enumerate(DetData.SignalArray):
		if (i == ind):
			continue
		if ( get_distance(s.pos, currSignal.pos) <= DetData.Var.neigh ):
			neighSignals.append(s)

	if( len(neighSignals) == 0 ):
		return False

	# Loop over those in limit
	for i,n in enumerate(neighSignals):
		if(n.snr_val > currSignal.snr_val):
			return True

	return False



# Extracts signal from map, and performs data analysis on it.
def process_signal(DetData, detectionMap, loop_counter):
	# Get index position of maximum SNR in snorm
	map_max = np.amax(detectionMap)
	map_size = detectionMap.shape[0]

	# Get coordinates of max in map
	index_max = get_max_ind(detectionMap)

	# Generate Mask
	rt = math.ceil(DetData.Var.res)
	polaire_outputs = idl.polaire2(rt, width=map_size,
						center_x=index_max[0], center_y=index_max[1],
						verbose = False)
	mask = polaire_outputs[2]

	# Calculate subwindow position indices in map
	subw_centre = (DetData.Var.subw_sz - 1)/2.0
	wx = [ int(index_max[0]-subw_centre) , int(index_max[0]+subw_centre) ]
	wy = [ int(index_max[1]-subw_centre) , int(index_max[1]+subw_centre) ]
	
	# Crop out signal subwindow on maps
	wmap = crop_window(DetData.Map.smap, mask, wx, wy)
	wnorm = crop_window(DetData.Map.snorm, mask, wx, wy)
	wdev = crop_window(DetData.Map.sdev, mask, wx, wy)
	wflux = None
	if(not DetData.Map.noflux):
		wflux = crop_window(DetData.Map.fmap, mask, wx, wy)

	
	# Create new signal class
	newSignal = SignalClass(wmap=wmap, wnorm=wnorm,wdev=wdev,wflux=wflux)
	newSignal.fit_gaussian(DetData, detectionMap)
	newSignal.fit_gaussian_flux(DetData)
	DetData.SignalArray.append(newSignal)

	# Save 3D plots
	if( DetData.Var.plot3d ):
		plot_3d_window_gaussian(sz = DetData.Var.subw_sz, 
			fwhm = [newSignal.fwhm_x, newSignal.fwhm_y],
			window = newSignal.maps.smap,
			path = DetData.Var.path,
			filename = "signal_"+str(loop_counter) )

	# ====== FLUX ======
	if(DetData.Map.noflux):
		return
	

	#---
	# Flux calculations
	# Apply Gaussian fit on flux map
	# Check if good fit and save flux value
	

# Crops out current signal from map for next detection
def crop_out_signal(DetData, detectionMap, loop_counter):
	rt = math.ceil(DetData.Var.res)
	img_sz = detectionMap.shape[0]
	imax = get_max_ind(detectionMap)
	# Get mask at signal position
	polaire_outputs = idl.polaire2(rt, width=img_sz, 
									center_x=imax[0], center_y=imax[1],
									verbose = False)
	mask_erase = polaire_outputs[2]
	# Multiply detection map by opposite of mask to remove signal
	detectionMap *= (1 - mask_erase)
	return detectionMap

# Translates numerical position flag to descriptive name
def pos_flag_name(flag):
	if(flag == -1):
		return "one-bad-one-good"
	elif(flag == 0):
		return "firm-detection"
	elif(flag == 1):
		return "both-bad-fit"
	elif(flag == 2):
		return "close-to-edge"
	elif(flag == 3):
		return "close-to-high-SNR-signal"
	elif(flag == 4):
		return "low-SNR"
	else:
		return "none"

# Saves detection data onto csv
def save_csv(DetData):
	# POSITIONS
	# Making arrays of class attributes
	coord_x = []
	coord_y = []
	snr = []
	snrerr = []
	flag = []
	flagname = []

	for s in DetData.SignalArray:
		coord_x.append(s.pos[0])
		coord_y.append(s.pos[1])
		snr.append(s.snr_val)
		snrerr.append(s.snr_err)
		flag.append(s.pos_flag)
		flagname.append(pos_flag_name(s.pos_flag))
		# FLAG APPEND

	# Save positions to Pandas dataframe
	d = {'x': coord_x, 'y': coord_y,
		'snr': snr, 'snr_err': snrerr,
		'flag': flag, 'flag_name': flagname}
	df = pd.DataFrame(data=d)
	df.to_csv(DetData.Var.path + "detections.csv")

	if( DetData.Var.verbose ):
		print("\n")
		print(" POSITION INFORMATION:")
		print(df)

	# FLUXES
	if(DetData.Map.noflux):
		return

	coord_x = []
	coord_y = []
	flux = []
	fluxerr = []
	contrast = []
	contrast_err = []
	fluxflag = []
	fluxflagname = []

	for s in DetData.SignalArray:
		coord_x.append(s.pos[0])
		coord_y.append(s.pos[1])
		flux.append(s.flux_val)
		fluxerr.append(s.flux_err)
		contrast.append(s.contrast)
		contrast_err.append(s.contrast_err)
		fluxflag.append(s.pos_flag)
		fluxflagname.append(pos_flag_name(s.pos_flag))
		# FLAG APPEND

	# Save positions to Pandas dataframe
	d = {'x': coord_x, 'y': coord_y,
		'flux': flux, 'flux_err': fluxerr,
		'contrast (mag)': contrast, 'contrast_err': contrast_err,
		'flag': fluxflag, 'flag_name': fluxflagname}
	df = pd.DataFrame(data=d)
	df.to_csv(DetData.Var.path + "fluxes.csv")

	if( DetData.Var.verbose ):
		print("\n")
		print(" FLUX INFORMATION:")
		print(df)


# Constructs map with detections and plots it
def save_detmap(DetData):
	plt.imshow(DetData.Map.snorm)
	plt.title('Detections in SNR Map')
	
	o_patch = mpatches.Patch( color='orange', label='Only one fit (SNR or Flux) is good')
	g_patch = mpatches.Patch( color='w', label='Firm detection')
	k_patch = mpatches.Patch( color='c', label='Close to high SNR signal')
	y_patch = mpatches.Patch( color='y', label='Close to edge')
	r_patch = mpatches.Patch( color='r', label='Bad SNR & Flux Fits')
	w_patch = mpatches.Patch( color='k', label='Low SNR')

	plt.legend(handles=[g_patch, o_patch, r_patch, k_patch, y_patch, w_patch], 
				loc = 'upper center',
				bbox_to_anchor=(1.05, 1),
				prop = {'size': 8},
				ncol = 1, fancybox = True, shadow = True)

	for s in DetData.SignalArray:
		color = 'r'
		if s.pos_flag == -1:
			color = 'orange'
		# Good fit
		elif s.pos_flag == 0:
			color = 'w'
		# Bad fit
		elif s.pos_flag == 1:
			color = 'r'
		# Too close to edge
		elif s.pos_flag == 2:
			color = 'y'
		# Too close to other signal
		elif s.pos_flag == 4:
			color = 'c'
		# Too close to high SNR neighbour
		elif s.pos_flag == 3:
			color = 'k'
		# Print circle
		circ = plt.Circle((s.pos[0], s.pos[1]),10, fill=False, color=color)
		plt.gcf().gca().add_artist(circ)

	plt.savefig(DetData.Var.path + "detections.png")
	plt.show()


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

	smap		: 2D float array
				SNR-map provided by ANDROMEDA
				-> For detection and subpx astrometry.

	snorm 		: 2D float array, optional.
				Normalised SNR map provided by ANDROMEDA

	snrdev 		: 2D float array, optional.
				Standard Deviation SNR map.

	fmap 		:  2D float array  
				FLUX-map provided by ANDROMEDA.
				-> For photometry.

	fdev 		: 2D float array
				Map of the standard deviation  
				of the estimated flux provided by ANDROMEDA.
				-> For error on contrast estimation
				-> And for detection limit computation.

	fndev 		:2D float array
				Map of normalised standard deviation of flux.

	limit		: integer (default = None), optional.
				threshold chosen so that detections are signal.
				found in the SNR-map that are above this threshold.

	subw_sz		: even integer [pixels]
				size of the 'subimages' fully containing the signal
				and in which the 2D-Gaussian fit are performed.

	pxscale		: float [arcsec/pixelSize]
				equivalent size of one pixel in arcsec 
				(given by the camera used for imagery).

	oversampl	: float [pixels]
				oversampling of the images given by the camera.

	neigh		: integer [pixels], optional
				maximum distance between the main lobe and the
				tertiary (positive) lobe in the pattern 
				of a planet signature made in the SNR-map. 
				If the main lobe has a very high SNR, it 
				may happen that the tertiary lobe is above 
				the threshold and thus detected as a companion
				whilst being an artefact.

	owa			: float [lambda/D], optional
				External limit of the data in the images (radius)

	iwa			: integer [lambda/D], optional
				Internal limit of the data in the images (radius),
				from which planetary signals are sought.

	telescope_type	: integer, optional
					Defines constraints for Gaussian Fit:
					Available: naco, nirc2, stim, other

	save_plots	: bool, optional
				If true, generates 3d plots of signals and their gaussian fits
				on 'detections' folder.

	verbose		: bool (default = True), optional.
				Displays coments while compiling the detection.

	max_signals : int (default = 10), optional
				Maximum number of detections to process.

	path 		: string (default = "./", current directory), optional
				Path to star data, where output files will be generated


	IDL Functions Required by original IDL script:
		im_contok.pro
		mpfit2dpeak.pro (and all subfunctions needed)
		gauss_calc.pro (cantallf)
		saveimage.pro 
		detection_limit.pro (cantallf)
		tvwin.pro (ONERA);   
	

	"""
def detection_andromeda(
			smap = None, snorm = None, sdev = None,
			fmap = None, fdev = None, fndev = None,
			limit = 0, subw_sz = 0, pxscale = 0, oversamp = 0,
			neigh = 0, owa = 0, iwa = 0,
			verbose = True, save_plots = False,
			telescope_type = "nirc2", max_signals = 10, path = "./"):

	# Initialising classes
	Maps = MapsClass(smap=smap, snorm=snorm, sdev=sdev,
					fmap=fmap, fdev=fdev, fndev=fndev,
					verbose=verbose)
	Vars = VarsClass(limit=limit, pxscale=pxscale, oversamp=oversamp,
					owa=owa, iwa=iwa, subw_sz=subw_sz,
					map_sz=smap.shape[0], neigh=neigh,
					verbose = verbose, plot3d = save_plots, path=path)
	GaussParams = GaussParamsClass(Vars.res, Vars.subw_sz, telescope_type)
	DetData = DetDataClass(Maps, Vars, GaussParams)

	# TO-DO
	# Make copy of snorm to edit and crop looking for signals
	detectionMap = np.copy(DetData.Map.snorm)

	# Check that at least a value is above threshold in snorm
	if( np.amax(detectionMap) < DetData.Var.limit ):
		print(" Warning: No signals found above threshold! Try a lower threshold.")
		return

	# LOOP FOR SIGNALS (or before reaching PMAX):
	for r in range(0, max_signals):
		
		# Check if all signals above threshold have already been found and cropped out
		if( np.amax(detectionMap) < DetData.Var.limit ):
			if(verbose):
				print(" Total signals found:", r)
			break

		if(verbose):
			print("\n Signal", r)
		process_signal(DetData, detectionMap, r)
		crop_out_signal(DetData, detectionMap, r)

	# Set positional flags for signals
	for i,s in enumerate(DetData.SignalArray):
		s.set_pos_flag(DetData, i)

	# Set flux flags

	# Save results in csv files
	save_csv(DetData)
	# Build detections map and plot
	save_detmap(DetData)


# =================================================================
# 		BEGINNING - MAIN
# =================================================================

def parse_args():

	parser = argparse.ArgumentParser()

	parser.add_argument("-v", "--verbose", action='store_true',
		help="Displays extra information during companion detection. Default = False")

	parser.add_argument("-p", "--plot3d", action='store_true',
		help="Generates 3D plots of signal and gaussian fits. Default = False")

	parser.add_argument("-g", "--gauss", type=str,
		help=("Specifies set of gaussian constraints to use on signals."
				" Names: naco, nirc2, stim, other"))

	parser.add_argument("-m","--maxdet", type=int,
		help="Max limit of detections in the map. Default = 10")


	args = parser.parse_args()

	# ---Flag for constraints
	# naco 	:NaCo-VLT
	# nirc2 :NIRC2-Keck
	# stim 	:STIM on NIRC2-Keck

	if( not args.gauss ):
		args.gauss = "nirc2"

	if( args.maxdet and args.maxdet <= 0):
		print(" Error: --maxdet must be an integer greater than zero")
		exit(-1)

	if( args.maxdet == None):
		args.maxdet = 10

	return args.verbose, args.plot3d, args.gauss, args.maxdet




# This is where the fun begins
def main():
	print("\n === Andromeda Detection Algorithm ===")

	# COMMAND LINE ARGUMENTS
	verbose, plot3d, obs_type, maxdet = parse_args()

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
		print("\n")
		print(" Input folder name where data is located: ")
		print(" e.g. 	'datasets/hr_8799/' or '.' for current folder.")
		path = input(" > ")
		print("\n")

		# Support for Windows + Unix systems: use / always for paths
		path.replace("\\","/")
		if(path[-1] != '/'):
			path += '/'
		
		print(" Input star name: ")
		print(" e.g. 	ANDR_snrmap_HIP_4067.fits --> HIP_4067")
		starname = input(" > ")
		print("\n")

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
	snr_stdev = None
	flux_map = None
	stddev_flux_map = None
	stddev_flux_norm = None

	if (option == 1):
		snr_stdev = path + "andr_snrmap_stddev_" + starname + ".fits"
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
		snr_norm = None
		flux_map = None
		stddev_flux_map = None
		stddev_flux_norm = None
		# Generate path
		slash = starpath.rfind("/")
		if(slash == -1):
			path = "./"
		else:
			path = starpath[0:slash+1]

	# OVERSAMPLING CALCULATION
	# Pixscale for K-band
	Kp_wavelength = 2.124E-6
	pixscale_keck = 0.00953
	diam_tel = 10
	pixscale = pixscale_keck * 1E3
	pixscale_nyq = 1/2 * (180*3600*1E3 / np.pi) * (Kp_wavelength / diam_tel)
	oversampling = pixscale_nyq / pixscale

	# THRESHOLD
	# For Naco/VLT ~= 5
	# For Nirc2/Keck ~= 10
	# For STIM ~= 1 (varies depending of noise)
	threshold = 0.0
	while(True):
		th = input(" Input threshold: ")
		try:
			th = float(th)
		except ValueError:
			print(" Error: threshold must be a float")
			continue

		if (th <= 0):
			print(" Error: threshold must be greater than zero")
			continue
		threshold = th
		break

	print(" Running...")

	detection_andromeda(smap=snr_map, snorm=snr_norm, sdev = None,
			fmap = flux_map, fdev = stddev_flux_map, fndev = stddev_flux_norm,
			limit = threshold, subw_sz = 0,
			pxscale = pixscale, oversamp = oversampling,
			neigh = 0, owa = 0, iwa = 0,
			verbose = verbose, telescope_type = obs_type, save_plots = plot3d,
			max_signals = maxdet, path = path)


if __name__ == '__main__':
	main()
	print("End of program")
	exit(0)