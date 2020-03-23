Stellar Ages Code README
Code courtesy of Trevor David (tdavid@caltech.edu)

RUNNING THE CODE:
1) Open terminal and navigate to the directory that has the required files inside. 
2) Type 'jupyter notebook' into the terminal and Jupyter tab should open in Safari.
3) Open 'CAMD_age_analysis.ipynb' in Jupyter
4) Select 'Cell' in the taskbar and then click 'Run All' to run the code, or press the 'Run' button to run the selected block of code individually. Once the code is running, it will print the progress at the bottom of the page.
5) The CMDs and their corresponding age histograms are outputted into the 'figures' folder. The solid red line is the median and the dotted lines are the 16th and 84th percentiles.
6) The ages with their corresponding errors are outputted as a tex file 'ages.tex'.
7) The raw histogram data are outputted to the 'agedist' folder.

GETTING DATA FROM VizieR
1) Go to http://vizier.u-strasbg.fr/viz-bin/VizieR. Search for Gaia in the catalog search. Select Gaia DR2 from the results and click 'Query selected catalogs'. Select the top option on the list 'I/345/gaia2'
3) Now you are ready to input your targets and retrieve Gaia data. You can either search targets individually or as a batch.
4) To batch search select 'List of targets' tab at the top of the page. You can input your list as a simple text file, 'targets.txt' is an example. 
5) Tick the 'add your input as first column' box. Change target dimension to 5 arcsec. Change the table output to ';-Seperated-Values' in the Preferences sidebar on the left of the page.
6) Select the following fields from the simple constraint list:
	RA_ICRS	
	e_RA_ICRS
	DE_ICRS	
	e_DE_ICRS
	Source	
	Plx	
	e_Plx	
	pmRA	
	e_pmRA	
	pmDE	
	e_pmDE	
	gofAL	
	chi2AL	
	epsi	
	sepsi	
	Dup	
	Gmag	
	e_Gmag	
	BPmag	
	e_BPmag	
	RPmag	
	e_RPmag	
	BP-RP	
	RV	
	e_RV	
	Teff	
	AG	
	E(BP-RP)
	Rad	
	Lum	
4) Once you have uploaded your target list and selected the constraints you want, click 'Submit' at the top of the page. The table should open in a new window. Copy the whole page and save it as a .tsv file in to the same folder as the code.
5) The data is now ready to be used in the code.
