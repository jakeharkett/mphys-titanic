# Age Analysis Code

This program calculates the age of a star with upper and lower errors, based on its Gaia magnitudes and parallax, and the MIST Isochrone models.

## Inputs

This program requires the following files:

### MIST models
MIST - MESA Isochrones and Stellar Tracks. Download from http://waps.cfa.harvard.edu/MIST/model_grids.html. Isochrones: Synthetic Photometry, v/vcrit=0.4, UBV(RI)c + 2MASS JHKs + Kepler + Hipparcos + Tycho + Gaia (116MB). Choose the option with zero metallicity " MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_basic.iso"

### Star Gaia data
Retrieve from VizieR: https://vizier.u-strasbg.fr/viz-bin/VizieR. Choose 'Gaia' mission, and download the fields Plx, e_Plx, Gmag, e_Gmag, BPmag, e_BPmag, RP_mag, and e_RPmag.

## Outputs
Saves the ages in both a text file 'ages.txt', and a tex file formatted as a LaTeX table 'ages.tex'. It also saves histogram data in the 'agedits' folder, and both histograms and CMDs for each star in the 'figures' folder.

## Possible Errors
* Error: file -- not found
    Star data file not found given input filename.

* Error: star data file has unexpected formatting.
	Expected column names not found in star data file.
	Please, check requirements above for star data formatting.
	Use a semicolon delimiter, with the first row being the column names,
	and rows below the corresponding data for each star.

* Error: models file 'MIST_Gaia_vvcrit0.4.iso.cmd' not found.
	The file with the MIST CMD models was not found.
	Please, ensure that is named 'MIST_Gaia_vvcrit0.4.iso.cmd'
	and is placed in same directory as this script.

* Error on star --: '--' field is not a number.
	One of the star's data fields could not be interpreted as a float.
	Either because it's missing or has an invalid value.

* Error on star --: parallax distribution has negative values (errors too large).
	Parallax of star has uncertantities that are too large,
	causing the generated distribution to have values lower or equal to zero.
	Try again with a smaller number of points in Monte Carlo simulation.
	
* Error on star ---: could not interpolate to model.
	Could not interpolate star data with MIST models for an unknown reason.
	Most likely, the generated distributions overran the model bounds.
	Try again with a smaller number of points in Monte Carlo simulation.
