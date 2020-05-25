# Detection Algorithm

This algorithm finds signals above a threshold in a given image and fits gaussian curves onto them. Place the python script anywhere, and run through command line interface. Then, choose whether to analyse the Andromeda Algorithm outputs, or just a single image.

## Input
To analyse Andromeda algorithm outputs, first input the path to the outputs from the location of the script (if it's in the current folder, simply input '.'). Then input the name of the star as written in the image filenames (e.g. 'hr_8799'). For analysing a single image, simply provide the location of said image. Finally, input the analysis threshold, which varies depending on the image to analyse. On the VLT SNR map given by Andromeda, the default is 5. On the Keck SNR map, the default is 10. On STIM maps, this value can be as low as 0.5.

## Output
It outputs an image with circles around the detections, colour coded for good or bad fits, low SNR, or proximity to other signals or the edges of the image. This image must be saved manually when generated. A 'detections.csv' file is generated, listing all the positions, SNR, and flags of every detection; and another 'fluxes.csv' file is generated listing the fluxes and uncertainties of the detections.

Furthermore, when setting the parameter "save_plots" to true in the code, every detected signal is plotted along with its gaussian fit and saved into a folder 'detections'.
