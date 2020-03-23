# Detection Algorithm

The algorithm must be called from the python file.
Place detection maps on same directory.
Ideally outputs of Andromeda Post-processing algorithm, but any map will work with an approrpiate
threshold.

It outputs signals found in a file 'detections.csv', along with its fluxes (if flux maps are
provided) in another file 'fluxes.csv'.
Optionally, it can also output 3D plots of the gaussian fits compared to the signa¨els, but it
requires an extra folder 'detections'.
