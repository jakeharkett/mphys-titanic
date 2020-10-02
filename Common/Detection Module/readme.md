# Detection Algorithm

This algorithm finds signals above a threshold in a given image and fits gaussian curves onto them. Place the python script anywhere, and run through command line interface. Then, choose whether to analyse the Andromeda Algorithm outputs, or just a single image.

Developed by Faustine Cantalloube in IDL. Translated to Python by Jorge Fernandez (jorgefz.fernandez@gmail.com)

Implemented from version:
Revision: 1.5, Date: 2018/06/25 

Check the paper Cantalloube et al. (2015) for more info: https://www.aanda.org/articles/aa/abs/2015/10/aa25571-14/aa25571-14.html

## Running the script

For default options, simply run the script in Python with one of the following commands:

`		py detection_andromeda.py`

`		python3 detection_andromeda.py`

## Inputs and Outputs
Once the scipt is running, you will be prompted to schoose between analysing an Andromeda Output set (1) or a single image (2).

Choosing to analyse the Andromeda Output, you will be askedto input the folder in which the images can be found. If the images are in the Then, input the starname as seen in the files' names: e.g.	ANDR_snrmap_HR_8799.fits has star name 'HR_8799'.

Choosing to analyse a single image will simply prompt you to type the path to the image.

Finally, you will be asked for a threshold. Only signals with a Signal-to-Noise ratio larger than this number will be sought in the Normalized SNR map (or the regular SNR map if the latter is missing).

Running the script will generate a map 'detections.png' wherein the signals are circled and colour-coded, as well as a file 'detections.csv' with signal data and 'fluxes.csv' with their flux data (if flux map was provided). These files will be generated on the path provided, where the input dataset or image is located.


## Command Line Arguments

The script may be called with a number of command line arguments. To check the available commands, run either of the following:
		
```
py detection_andromeda.py -h
py detection_andromeda.py --help
```

The commands are the following:

```
		usage: newdet.py [-h] [-v] [-p] [-g GAUSS] [-m MAXDET]

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
```

Examples:

To run the script with verbose messages (extra info):
```
			py detection_andromeda.py -v
			py detection_andromeda.py --verbose
```

To run the script with the verbose flag and save 3D plots:
```
			py detection_andromeda.py -v -p
			py detection_andromeda.py --verbose --plot3d
```

To run the script with a constraint set 'naco':
```
			py detection_andromeda.py -g naco
			py detection_andromeda.py --gauss naco
```

To run the script with 50 max detections:
```
			py detection_andromeda.py -m 50
			py detection_andromeda.py --maxdet 50
```
