**************************************NOTES****************************************

PROJECT SCRIPT 12

Script last updated 21/01/2020
Works for VIP 0.9.11 in Python3 environment (current version on analysis)

For any help, contact:
	Jorge Fernández
	MPhys student with Sasha, 2019/2020
	jorgefz.fernandez@gmail.com

**************************************INSTRUCTIONS**********************************

Place script on the folder where the science images are.
Run with python.
Enter the name of the star to run the script, e.g. HIP_99542. 


Ignore this for the moment:
	--Please upload contrcurve.py, cosmetics.py, fakecomp.py and shapes.py when running this script (avaliable in DangerZone)



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


----------------------------------------

DEVELOPMENT CONTINUED BY JORGE FERNANDEZ

Contact: jorgefz.fernandez@gmail.com


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

09/10/19
	Fixed some Contrast curves problems 
		Still needs optimal number of PCs
	PCA function with full output doesn't have the expected number of outputs
		-Need residual cube for STIM maps
		
15/11/19
	Changed name to Projectscript12
	Removed nested loops that would try different values to recenter star
	Added option to choose between using 2D Gaussian fit or manual fit when recentering
	
18/11/19
	Implementing Andromeda

19/11/19
	Reordered post-processing algorithms into functions
	Obtained optimal PCs from FFPCA output
	Implemented Annular PCA
	Implemented STIM map using residual cube from Annular PCA
	Implemented STIm map using residual cube from LLSG

20/01/2020
	Updated Andromeda Oversampling Factor

21/01/2020
	Fixed Contrast Curves with James Davies' work.

23/03/2020
	Reworked Recentering function
	Added more recentering algorithms

Script designed to build a data cube from fits files with the aid of the VIP module
in python, reduce the images using the PCA and LLSG method of imaging reduction. 
Script also then calculates the flux from the exposure times of the fits files and 
then generates contrast curves based on the data cube produced and this value of flux.


STUCTURE

	1. Packages
	
	2. Global variables
	
	3. Functions

	4. Main Functions

	5. Image Processing

