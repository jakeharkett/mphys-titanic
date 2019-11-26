"""
            RDI.py
            Author: Brendan Wilby
            Created 26/11/2019

            A script to perform ADI and/or RDI on a science cube.
            If a reference cube is provided, this script will perform RDI. If not it will default to ADI.

            USAGE:
                python RDI.py [science cube] [refcube=None]

            TODO:
                - Does (my) analysis need upgrading to python 3 + vip.0.9.11?????
                - How to input the science cube?
"""


import vip_hci as vip
from astropy.io import fits
import numpy as np
import sys

HEADER_NAME_PARANG = "PARANG"
PSF_SOURCE = "psf.fits"
PSF_XY = [95, 111]


"""
=================================================================================================
                    HELPER FUNCTIONS
=================================================================================================
"""

def ReadAngles(target_name, save=False):
    file_names_path = target_name + "_filenames.txt"
    try:
        file_names = np.loadtxt(file_names_path, dtype=str)
    except NameError:
        print("Could not load file %s" %file_names_path)
        sys.exit()

    num_images = len(file_names)
    angles = np.zeros(num_images)

    for i in range(num_images):
        hdulist = fits.open(file_names[i], ignore_missing_end=True)
        angles[i] = hdulist[0].header[HEADER_NAME_PARANG]

    if save:
        save_filename = target_name + "_angles.txt"
        np.savetxt(save_filename, angles)
        print("========================================================================================================")
        print("Saved angles to file %s" %save_filename)
        print("========================================================================================================")

    print("========================================================================================================")
    print("Successfully loaded %d angles for target %s" %(num_images, target_name))
    print("========================================================================================================")

    return angles


"""
=================================================================================================
                    MAIN FUNCTION
=================================================================================================
"""
if __name__ == "__main__":
    num_args = len(sys.argv)
    if num_args < 2 or num_args > 3:
        print("Invalid command entered. Usage: python RDI.py [science cube] [refcube=None]")
        sys.exit()

    ref_cube = None
    science_cube = None
    science_cube_path = ""
    ref_cube_path = None

    # Load target
    target_name = sys.argv[1]

    try:
        science_cube_path = "centeredcube_" + target_name + ".fits"
        hdulist = fits.open(science_cube_path, ignore_missing_end=True)
        science_cube = hdulist[0].data
        print("========================================================================================================")
        print("Loaded science cube %s of shape: (%d, %d, %d)" %(science_cube_path, science_cube.shape[0], science_cube.shape[1], science_cube.shape[2]))
        print("========================================================================================================")
    except NameError:
        print("Could not load file %s" %science_cube_path)
        sys.exit()

    # Load in the reference cube if command line argument [refcube] is specified
    if num_args == 3:
        ref_cube_path = sys.argv[2]

        try:
            hdulist = fits.open(ref_cube_path, ignore_missing_end=True)
            ref_cube = hdulist[0].data
            print("========================================================================================================")
            print("Loaded reference cube %s of shape: (%d, %d, %d)" %(ref_cube_path, ref_cube.shape[0], ref_cube.shape[1], ref_cube.shape[2]))
            print("========================================================================================================")
        except NameError:
            print("File %s does not exist." %ref_cube_path)
            sys.exit()

    angle_list = ReadAngles(target_name, True)

    """
        WHAT ARE THESE??????????????????????
    """
    inner_rad_rdi = 0.1         # Guessed value
    PCA_comps = 1               # how to find this?

    try:
        psf = vip.fits.open_fits(PSF_SOURCE, n=0, header=False, ignore_missing_end=True, verbose=False)
    except NameError:
        print("Could not load file at %s" %PSF_SOURCE)
        sys.exit()

    gauss=vip.var.fit_2dgaussian(psf, crop=True, cropsize=30, cent=(PSF_XY[0], PSF_XY[1]), full_output=False, debug=False)

    print(gauss[0:1])
    fwhm_x = gauss[0]
    fwhm_y = gauss[1]
    fwhm = np.mean([fwhm_x, fwhm_y])
    print(fwhm)

    mask_center_pixels = inner_rad_rdi * fwhm

    output_cube = vip.pca.pca_fullfr.pca(
        cube=science_cube, angle_list=angle_list, cube_ref=ref_cube,
        ncomp=1, svd_mode="lapack", scaling="spat-mean",
        mask_center_px=mask_center_pixels, source_xy=(PSF_SOURCE[0], PSF_SOURCE[1]),
        fwhm=fwhm, full_output=False, verbose=True)

    output_filename = target_name + "_RDI.fits"
    hdu_new = fits.PrimaryHDU(output_cube)
    hdu_new.writeto(output_filename, overwrite=True)
    print("Finished writing to output/%s." %output_filename)