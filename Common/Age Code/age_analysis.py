






"""
===============================================
	Colour Magnitude Diagram Age Analysis
===============================================

Original Author:
	Trevor David (tdavid@caltech.edu)

Further developed by Jorge Fernandez
	(jorgefz.fernandez@gmail.com)


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
		File that contains isochrones to plot
		Various ages and for various filters (Gaia mags, 2MASS mags...)
		Here we use Gaia mags, and ages 1, 10, 100 Myr, and 1 10 Gyr tracks
		Source:
			MIST - MESA Isochrones and Stellar Tracks
			https://arxiv.org/abs/1604.08592

	OUTPUTS:
		- Calculated ages and uncertainties for each star
		in files 'ages.tex' and 'ages.txt'
		- CMD plot with all input stars in 'global_cmd.png'
		- Individual CMD plots and age histograms in 'figures' folder
		- Histogram data in 'agedist' folder
	

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
	_mu   = 5.0*np.log10(_dist)-5.0           # Distance modulus
	_MG   = gmag-_mu
	return _MG

#def BPRP(bpmag, e_bpmag, rpmag, e_rpmag):
    #Again, assuming BP and RP errors are uncorrelated
#    _bp = np.random.normal(bpmag, e_bpmag, npts)
#    _rp = np.random.normal(rpmag, e_rpmag, npts)
#    return _bp-_rp







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
	    
	if len(xdata)>0:
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


def calculate_ages(mist, dr2):

	npts = int(1e6) #number of points to use in Monte Carlo simulation

	xmodel = mist['Gaia_BP_DR2Rev']-mist['Gaia_RP_DR2Rev']
	ymodel = mist['Gaia_G_DR2Rev']
	zmodel = 10.**mist['log10_isochrone_age_yr']/1.0e6

	MG = absolute_Gmag(dr2['Gmag'], dr2['Plx'])
	BPRP = dr2['BPmag']-dr2['RPmag']

	#Write results to TeX file in LaTeX format
	save_tex = open('ages.tex', 'w')

	#Write results to normal text file, tab-separated values for
	# name age +error -error
	save_txt = open('ages.txt', 'w')

	for i in range(len(dr2)):

		# Monte Carlo simulations to calculate age and uncertainties,
		# accounting for errors in photometry and parallax

		#Turning all data into gaussian distributions
		_gmag = np.random.normal(dr2['Gmag'][i], dr2['e_Gmag'][i], npts)
		_plx  = np.random.normal(dr2['Plx'][i], dr2['e_Plx'][i], npts)
		_MG   = absolute_Gmag(_gmag, _plx)
		_BP   = np.random.normal(dr2['BPmag'][i], dr2['e_BPmag'][i], npts)
		_RP   = np.random.normal(dr2['RPmag'][i], dr2['e_RPmag'][i], npts)
		_BPRP = _BP-_RP

		if np.isnan(_BPRP).any() or np.isnan(_MG).any():
			print(i, "Error on star", dr2['_1'][i], "(BPRP or MG)")
			save_tex.write("%s\n" % dr2['_1'][i])
			save_txt.write("%s\n" % dr2['_1'][i])
			continue

		# np.random.normal generates a normal distribution based on a median and two errors
		# with a number of points equal to 'npts', here set to 1e6

		# Interpolate gaussian points with isochrones to get age gaussian
		_agemyr = model_2dinterp(_BPRP, _MG, xmodel, ymodel, zmodel)

		if np.isnan(_agemyr).any():
			printf(i, "Error on star", dr2['_1'][i], "(AGEMYR)")
			save_tex.write("%s\n" % dr2['_1'][i])
			save_txt.write("%s\n" % dr2['_1'][i])
			continue

		# Runs through age distribution values looking for NaN, meaning interpolation
		# could not be completed
		    
		# Get median age and errors
		_agemed  = np.nanpercentile(_agemyr, 50.)
		_agelerr = np.nanpercentile(_agemyr, 16.)-_agemed
		_ageuerr = np.nanpercentile(_agemyr, 84.)-_agemed

		# Write the individual age distributions to file
		f = open('./agedist/'+dr2['_1'][i]+'-agedist.txt', 'w')
		for j in range(len(_agemyr)):
			f.write('%.1f\n' % (_agemyr[j]))
		f.close()

		# Write LaTeX table
		save_tex.write('%s\n' % (dr2['_1'][i]+' & '+str("{:.2f}".format(_agemed))+'^{'+str("{:+.2f}".format(_ageuerr))+'}_{'+str("{:.2f}".format(_agelerr))+'} \\\\'))
		save_txt.write('%s\n' % (dr2['_1'][i]+'\t'+str("{:.2f}".format(_agemed))+'\t'+str("{:+.2f}".format(_ageuerr))+'\t'+str("{:.2f}".format(_agelerr)) ) )

		# Make CMD for each individual source for spot-checking purposes
		camd_age_plot(mist, xdata=[np.nanmedian(_BPRP)], ydata=[np.nanmedian(_MG)], xerr=[np.nanstd(_BPRP)], yerr=[np.nanstd(_MG)], age_array=_agemyr, filename='./figures/'+dr2['_1'][i]+'-camd.jpg')

		# Print progress
		if (i+1)%2==0:
			print("{:.1f}".format(100.*(i+1)/len(dr2)),' % complete')

	save_tex.close()
	save_txt.close()



def main():

	print("\n\n")
	print(" ====== CMD Age Analysis =======\n")
	print(" Input filename with star data (data delimiter is semicolon)")
	filename = input(" Filename = ")

	dr2 = np.genfromtxt(filename, 
						skip_header=0, 
						delimiter=';', 
						names=True, 
						dtype=None, 
						encoding='utf8', 
						autostrip=True)

	# Calculating absolute magnitudes, and colour index BPRP (blue - red magnitudes)
	MG = absolute_Gmag(dr2['Gmag'], dr2['Plx'])
	BPRP = dr2['BPmag']-dr2['RPmag']

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
	
	calculate_ages(mist, dr2)

	print("end")




if __name__ == "__main__":
	main()