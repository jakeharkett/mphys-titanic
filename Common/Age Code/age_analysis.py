






"""
===============================================
	Colour Magnitude Diagram Age Analysis
===============================================

Original Author:
	Trevor David (tdavid@caltech.edu)

Further developed by Jorge Fernández
	(jorgefz.fernandez@gmail.com)

	
	Developed in Python 3

	INPUTS:
		- Input stars
		List of stars with semicolon delimiter and the following data columns:
			_1		name
			Plx		parallax
			e_Plx	parallax error
			Gmag	Green Gaia mag
			e_Gmag	Green Gaia mag error
			BPmag	Blue Gaia mag
			e_BPmag	Blue Gaia mag error
			RPmag	Red Gaia mag
			e_RPmag	Red Gaia mag error

		- MIST Isochrones datafile
		Filename = "MIST_Gaia_vvcrit0.4.iso.cmd"
		File that contains isochrones to plot
		Various ages and for various filters (Gaia mags, 2MASS mags...)
		Here we use Gaia mags, and ages 1, 10, 100 Myr, and 1 10 Gyr tracks
		Source:
			MIST - MESA Isochrones and Stellar Tracks
			https://arxiv.org/abs/1604.08592
		Webpage:
			http://waps.cfa.harvard.edu/MIST/model_grids.html
			"Isochrones: Synthetic Photometry" -> "v/vcrit=0.4" -> "UBV(RI)c + 2MASS JHKs + Kepler + Hipparcos + Tycho + Gaia (116MB)"
		Download:
			http://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_vvcrit0.4_UBVRIplus.txz
			-> MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_basic.iso

		
	OUTPUTS:
		- Calculated ages and uncertainties for each star
		in files 'ages.tex' and 'ages.txt'
		- CMD plot with all input stars in 'global_cmd.png'
		- Individual CMD plots and age histograms in 'figures' folder
		- Histogram data in 'agedist' folder
	
	
	HOW TO USE:

	1) Place 'age_analysis.py' and 'MIST_Gaia_vvcrit0.4.iso.cmd'
	in same directory.

	2) Place file with stellar data in same directory.

	3) Execute python script using command line interface.

	4) Input filename of star data.

	5) Let the program run.


	POSSIBLE ERRORS:
	
	- Error: file -- not found
	Star data file not found given input filename.

	- Error: star data file has unexpected formatting.
	Expected column names not found in star data file.
	Please, check requirements above for star data formatting.
	Use a semicolon delimiter, with the first row being the column names,
	and rows below the corresponding data for each star.

	- Error: models file 'MIST_Gaia_vvcrit0.4.iso.cmd' not found.
	The file with the MIST CMD models was not found.
	Please, ensure that is named 'MIST_Gaia_vvcrit0.4.iso.cmd'
	and is placed in same directory as this script.

	- Error on star --: '--' field is not a number.
	One of the star's data fields could not be interpreted as a float.
	Either because it's missing or has an invalid value.

	- Error on star --: parallax distribution has negative values (errors too large).
	Parallax of star has uncertantities that are too large,
	causing the generated distribution to have values lower or equal to zero.
	Try again with a smaller number of points in Monte Carlo simulation.
	
	- Error on star ---: could not interpolate to model.
	Could not interpolate star data with MIST models for an unknwon reaso.
	Most likely, the generated distributions overrun the model bounds.
	Try again with a smaller number of points in Monte Carlo simulation.



"""




#%matplotlib inline

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from scipy.interpolate import griddata
import os


mpl.rcParams['font.size']=16
mpl.rcParams['figure.dpi']=300
mpl.rcParams['lines.markersize'] = 3
mpl.rcParams['savefig.bbox'] = 'tight'

# Notes on individual sources
# HIP 14809: Co-moving companion at 0.5"?
# HIP 16563: Co-moving companion at 0.2"?
# G 34-23: Target has a parallax of 87.6 +/- 8.5 mas, from 2016AJ....151..160F
# HIP 104430: Has Gaia photometry but parallax and proper motions come from Hipparcos! (van Leeuwen 2007)


# Function takes relative green magnitude and parallax
# Returns absolute green magnitude
def absolute_Gmag(gmag, plx):
	#NOTE: Computing distance simply as 1/plx. There may be offsets or biases that one needs to worry about,
	#particularly for the closest targets (see papers by Stassun & Torres, Bailer-Jones, etc.)
	#Also note, assuming Gmag and Plx errors are uncorrelated here.

	_dist = 1.0e3/plx
	_mu   = 5.0*np.log10(_dist) - 5.0
	_MG   = gmag - _mu
	return _MG







# Plots CMD with all targets
def camd_plot(mist, xdata=None, ydata=None, xerr=None, yerr=None):

	log10ages = np.arange(6,11,1)
	isoc_labels = ['1 Myr', '10 Myr', '100 Myr', '1 Gyr', '10 Gyr']
	xisoc = mist['Gaia_BP_DR2Rev']-mist['Gaia_RP_DR2Rev']
	yisoc = mist['Gaia_G_DR2Rev']

	for i in range(len(log10ages)):
	    arg = np.where((mist['log10_isochrone_age_yr']==log10ages[i])&(mist['EEP']<700))
	        # np.where gets data from MIST file that matches log10ages[i] age, and 'EEP' is less than 700 (resolution?) 
	    plt.plot(xisoc[arg], yisoc[arg], label=isoc_labels[i])
	    
	if xdata.size>0:
	    plt.plot(xdata, ydata, 'ko', mfc='None')
	    

	plt.xlabel(r'G$_\mathregular{BP}$-G$_\mathregular{RP}$ (mag)')
	plt.ylabel(r'M$_\mathregular{G}$ (mag)')
	plt.gca().invert_yaxis()
	plt.legend(prop={'size':10})
	plt.savefig('global_cmd.png')
	plt.show()
	return

# PLots all individual CMDs with age histograms
def camd_age_plot(mist, xdata=None, ydata=None, xerr=None, yerr=None, age_array=None, filename=None):

	log10ages = np.arange(6,11,1)
	isoc_labels = ['1 Myr', '10 Myr', '100 Myr', '1 Gyr', '10 Gyr']
	xisoc = mist['Gaia_BP_DR2Rev']-mist['Gaia_RP_DR2Rev']
	yisoc = mist['Gaia_G_DR2Rev']

	plt.subplot(121)

	for i in range(len(log10ages)):
	    arg = np.where((mist['log10_isochrone_age_yr']==log10ages[i])&(mist['EEP']<700))
	    plt.plot(xisoc[arg], yisoc[arg], label=isoc_labels[i])
	    
	if len(xdata)>0:
	    plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='o', mfc='None', color='k')
	    

	plt.xlabel(r'G$_\mathregular{BP}$-G$_\mathregular{RP}$ (mag)')
	plt.ylabel(r'M$_\mathregular{G}$ (mag)')
	plt.gca().invert_yaxis()
	plt.legend(prop={'size':10})

	plt.subplot(122)
	plt.hist(age_array, bins=100, histtype='step', color='k', density=True)
	plt.axvline(np.nanmedian(age_array), color='r')
	plt.axvline(np.nanpercentile(age_array, 16.), color='r', ls=':')
	plt.axvline(np.nanpercentile(age_array, 84.), color='r', ls=':')
	plt.xlabel('Age (Myr)')
	plt.ylabel('Probability density')

	plt.gcf().set_size_inches(12,5)
	plt.tight_layout()
	if filename != None:
	    plt.savefig(filename, dpi=100)
	plt.close()
	return


# Interpolates star age between the two closest isochrones
def model_2dinterp(xdata, ydata, xmodel, ymodel, zmodel, interp_method='linear'):
	# Interpolate over a 2D XY model grid (xmodel, ymodel) to determine corresponding Z values (zmodel)
	# for input XY data (xdata, ydata)
	# Options for interp_method are  'linear', 'cubic', 'nearest'
	points = xmodel, ymodel
	zdata = griddata(points, zmodel, (xdata,ydata), method=interp_method)
	return zdata

def check_star(star, i, star_num):

	cols = ('Plx', 'e_Plx', 'Gmag', 'e_Gmag',
			'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag')

	for c in cols:
		# Check if value can be interpreted as float
		try:
			float(star[c])
		except ValueError:
			print(" Error on star %d/%d (%s): '%s' field is not a number" % (i+1, star_num, star['_1'], c))
			return False
		# Check if value is NaN
		if (star[c] != star[c]):
			print(" Error on star %d/%d (%s): '%s' field is not a number" % (i+1, star_num, star['_1'], c))
			return False

	return True


def calculate_ages(mist, dr2):

	npts = int(1e6) #number of points to use in Monte Carlo simulation

	xmodel = mist['Gaia_BP_DR2Rev']-mist['Gaia_RP_DR2Rev']
	ymodel = mist['Gaia_G_DR2Rev']
	zmodel = 10.**mist['log10_isochrone_age_yr']/1.0e6

	MG = absolute_Gmag(dr2['Gmag'], dr2['Plx'])
	BPRP = dr2['BPmag'] - dr2['RPmag']

	# Arrays for ages
	ages = []
	uerr = []
	lerr = []

	for i in range(dr2.size):

		# Monte Carlo simulations to calculate age and uncertainties,
		# accounting for errors in photometry and parallax

		if not check_star(dr2[i], i, dr2.size):
			ages.append(0)
			uerr.append(0)
			lerr.append(0)
			continue
		
		_gmag = np.random.normal(float(dr2['Gmag'][i]), float(dr2['e_Gmag'][i]), npts)
		_plx  = np.random.normal(float(dr2['Plx'][i]), float(dr2['e_Plx'][i]), npts)

		if (_plx <= 0).any():
			print(" Error on star %d/%d (%s): parallax distribution has negative values (errors too large)" % (i+1,dr2.size,dr2['_1'][i]) )
			ages.append(0)
			uerr.append(0)
			lerr.append(0)
			continue

		_MG   = absolute_Gmag(_gmag, _plx)
		_BP   = np.random.normal(float(dr2['BPmag'][i]), float(dr2['e_BPmag'][i]), npts)
		_RP   = np.random.normal(float(dr2['RPmag'][i]), float(dr2['e_RPmag'][i]), npts)
		_BPRP = _BP - _RP

		# np.random.normal generates a normal distribution based on a median and two errors
		# with a number of points equal to 'npts', here set to 1e6

		# Interpolate gaussian points with isochrones to get age gaussian
		_agemyr = model_2dinterp(_BPRP, _MG, xmodel, ymodel, zmodel)

		if np.isnan(_agemyr).any():
			print(" Error on star %d/%d (%s): could not interpolate to model." % (i+1,dr2.size,dr2['_1'][i]) )
			ages.append(0)
			uerr.append(0)
			lerr.append(0)
			continue

		# Runs through age distribution values looking for NaN, meaning interpolation
		# could not be completed
		    
		# Get median age and errors
		_agemed  = np.nanpercentile(_agemyr, 50.)
		_agelerr = np.nanpercentile(_agemyr, 16.) - _agemed
		_ageuerr = np.nanpercentile(_agemyr, 84.) - _agemed

		ages.append(_agemed)
		uerr.append(_ageuerr)
		lerr.append(_agelerr)

		# Write age distribution data to file
		f = open('./agedist/'+dr2['_1'][i]+'-agedist.txt', 'w')
		for j in range(len(_agemyr)):
			f.write('%.1f\n' % (_agemyr[j]))
		f.close()

		# Make CMD for each individual source for spot-checking purposes
		camd_age_plot(mist, xdata=[np.nanmedian(_BPRP)], ydata=[np.nanmedian(_MG)], xerr=[np.nanstd(_BPRP)], yerr=[np.nanstd(_MG)], age_array=_agemyr, filename='./figures/'+dr2['_1'][i]+'-camd.jpg')

		# Print progress
		print(" Star %d/%d (%s) processed." % (i+1,len(dr2), dr2['_1'][i]))

	return ages, uerr, lerr


def check_columns(columns):
	expected_cols = ('_1', 'Plx', 'e_Plx', 'Gmag', 'e_Gmag',
					'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag')
	for n in expected_cols:
		if (n not in columns):
			print(" Error: star data file has unexpected formatting.")
			print(" 	-> '%s' data column missing" % n)
			print(" Check documentation for expected columns in data file.")
			exit()


def main():

	print("\n\n")
	print(" ====== CMD Age Analysis =======\n")
	print(" Input filename with star data:")
	filename = input(" Filename = ")

	if not os.path.exists(filename):
		print(" Error: file '%s' not found" % filename)
		return


	dr2 = np.genfromtxt(filename, 
						skip_header=0, 
						delimiter=';', 
						names=True, 
						dtype=None, 
						encoding='utf8', 
						autostrip=True)

	# Check correct number of column names in star data file.
	check_columns(dr2.dtype.names)

	#With only one entry in star data file, lenth of array will be zero,
	# which makes the code fail.
	# Solve by duplicating entry, and naming the duplicate 'ignore'
	if (dr2.size < 2):
		dr2 = np.append(dr2, dr2)
		dr2[1][0] = 'ignore'

	# Calculating absolute magnitudes, and colour index BPRP (blue - red magnitudes)
	MG = absolute_Gmag(dr2['Gmag'], dr2['Plx'])
	BPRP = dr2['BPmag']-dr2['RPmag']

	model_filename = "./MIST_Gaia_vvcrit0.4.iso.cmd"
	if not os.path.exists(model_filename):
		print("Error: models file '%s' not found." % model_filename)
		exit()

	# Extracting MIST models
	mist = np.genfromtxt('./MIST_Gaia_vvcrit0.4.iso.cmd', skip_header=12, names=True)

	# Making directories
	try:
		os.mkdir('agedist')
		os.mkdir('figures')
	except FileExistsError:
		pass

	# Plot CMD with all stars included
	camd_plot(mist, xdata=BPRP, ydata=MG)

	# Plot Inidividual CMDs
	#camd_age_plot(mist, xdata=BPRP, ydata=MG, xerr=0.1, yerr=0.1, age_array=np.random.normal(100,10,int(1e4)))
	ages, uerr, lerr = calculate_ages(mist, dr2)


	#Write results to TeX file in LaTeX format
	save_tex = open('ages.tex', 'w')
	#Write results to normal text file
	save_txt = open('ages.txt', 'w')

	save_tex.write("%% Star \t Age \t Upper Error \t Lower Error\n")
	save_txt.write("# Star \t Age \t Upper Error \t Lower Error\n")
	for i in range(len(ages)):
		save_tex.write("%s & %.2f ^{ +%.2f }_{ -%.2f } \\\\ \n" % (dr2['_1'][i], ages[i], uerr[i], lerr[i]))
		#save_tex.write('%s\n' % (dr2['_1'][i]+' & '+str("{:.2f}".format(ages[i]))+'^{'+str("{:+.2f}".format(uerr[i]))+'}_{'+str("{:.2f}".format(lerr[i]))+'} \\\\'))
		save_txt.write("%s \t %.2f \t +%.2f \t -%.2f \n" % (dr2['_1'][i], ages[i], uerr[i], lerr[i]))
		#save_txt.write('%s\n' % (dr2['_1'][i]+'\t'+str("{:.2f}".format(ages[i]))+'\t'+str("{:+.2f}".format(uerr[i]))+'\t'+str("{:.2f}".format(lerr[i])) ) )
	save_tex.close()
	save_txt.close()
	

	print("Done!")




if __name__ == "__main__":
	main()