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