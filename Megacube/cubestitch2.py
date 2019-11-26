"""
    cubestitch2.py - Written by Brendan Wilby November 2019

    USAGE:

    python cubestitch.py -create EPOCH_NAME
    python cubestitch.py -align EPOCH_NAME
    python cubestitch.py -angles EPOCH_NAME
"""


import sys
import math
import vip_hci as vip
from vip_hci.preproc import cube_recenter_2dfit, cube_recenter_via_speckles
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits

CUBES_FOLDER = "output/"
EPOCH_FOLDER = "targets/"
BAD_FRAMES_FOLDER = "badframes/"
CUBE_STITCHDATA_FOLDER = "stitch_data/"
TRIMMED_SUFFIX = "_trimmed"
CUBE_SUFFIX = "_megacube"
ALIGNED_SUFFIX = "_aligned"
PSF_SOURCE = "psf.fits"
PSF_XY = [95, 111]
IMG_DIMS = 1023


"""
    ============================================================================
                                HELPER FUNCTIONS
    ============================================================================
"""

def print_usage():
    print("USAGE:")
    print(" -create EPOCH_NAME  Create new cube from folder EPOCH_NAME")
    print(" -align  CUBE_NAME   Align cube at path CUBE_NAME")
    print(" -angles EPOCH_NAME  Extract parallactic angles from EPOCH_NAME .fits headers")

def resize_img(img):
    return img[:, :IMG_DIMS, :IMG_DIMS]

def get_centre_cubes(epoch, targets):
    return [epoch + "/" + target.replace("_", "") + "/centeredcube_" + target + ".fits" for target in targets]

def count_images(epoch):
    img_count = 0
    targets = load_targets(epoch)
    paths = get_centre_cubes(epoch, targets)

    for path in paths:
        try:
            hdulist = fits.open(path)
            hdu = hdulist[0]
            img_count += hdu.data.shape[0]
        except NameError:
            print("Path %s not found" %path)

    return img_count

def load_targets(epoch):
    targets = []
    path = EPOCH_FOLDER + epoch + ".txt"

    try:
        targets = np.loadtxt(path, dtype="S")
    except NameError:
        print("Could not find the path specified at " + path)

    return targets

def save_cube(cube, epoch, prefix):
    file_name = epoch + prefix + ".fits"

    hdu_new = fits.PrimaryHDU(cube)
    hdu_new.writeto("output/" + file_name, overwrite=True)
    print("Finished writing to output/%s." %file_name)

"""
    ============================================================================
                                COMMAND FUNCTIONS
    ============================================================================
"""

def create_cube(epoch_name):
    try:
        cube_data = open(CUBE_STITCHDATA_FOLDER + epoch_name + "stitch_data.csv", "w")
        print("Creating cube from epoch " + epoch_name)
        img_count = 0

        num_images = count_images(epoch_name)
        megacube = np.zeros((num_images, IMG_DIMS, IMG_DIMS))
        print("Created blank cube for %s of size %d images" %(epoch_name, num_images))

        targets = load_targets(epoch_name)
        paths = get_centre_cubes(epoch_name, targets)

        print("Beginning stitching of cubes...")

        for index, path in enumerate(paths):
            hdulist = fits.open(path)
            hdu = hdulist[0]
            num_images = hdu.data.shape[0]

            print("Loaded " + path)

            if hdu.data.shape[1] != IMG_DIMS:
                hdu.data = resize_img(hdu.data)
                print("Resized %s to shape (%d, %d)" %(path, IMG_DIMS, IMG_DIMS))

            for i in range(img_count, (img_count + num_images)):
                megacube[i, :, :] = hdu.data[i - img_count, :, :]
            
            img_count += num_images
            cube_data.write(targets[index] + ", " + str(num_images) + ", " + str(img_count) + "\n")


        print("Finished stitching")
        save_cube(megacube, epoch_name, CUBE_SUFFIX)
        cube_data.close()
    except NameError:
        print("Error: epoch %s not found" %epoch_name)  

def align_cube(cube_path_name):
    path = CUBES_FOLDER + cube_path_name + ".fits"
    print("Aligning cube")

    psf = vip.fits.open_fits(PSF_SOURCE, n=0, header=False, ignore_missing_end=True, verbose=True)
    gauss = vip.var.fit_2dgaussian(psf, crop=True, cropsize=30, cent=(PSF_XY[0], PSF_XY[1]), debug=True, full_output=False)

    fwhm_x = gauss[0]
    fwhm_y = gauss[1]
    fwhm = np.mean([fwhm_x, fwhm_y])
    fwhm_int = int(math.ceil(fwhm))

    hdulist = fits.open(path)
    hdu = hdulist[0]
    cube = hdu.data

    cube = cube_recenter_2dfit(cube, (512, 512), fwhm_int, subi_size=fwhm_int, model="gauss", nproc=1, debug=False, negative=True, full_output=False, plot=False)

    save_cube(cube, cube_path_name, ALIGNED_SUFFIX)
    


def remove_frames(epoch_name):
    frames_file = BAD_FRAMES_FOLDER + epoch_name + ".txt"
    bad_frames = []
    cube_path = CUBES_FOLDER + epoch_name + CUBE_SUFFIX + ".fits"

    try:
        bad_frames = np.loadtxt(frames_file, dtype=int)

    except NameError:
        print("Could not load file %s" %frames_file)
    except TypeError:
        print("No frames found in badframes file: %s" %frames_file)

    num_bad_frames = len(bad_frames)

    try:
        hdulist = fits.open(cube_path)
        hdu = hdulist[0]
        num_images = hdu.data.shape[0]

        print("Num images: %d" %num_images)
        print("Bad images: %d" %num_bad_frames)
        print("New images: %d" %(num_images - num_bad_frames))

        new_cube = np.zeros((num_images - num_bad_frames, IMG_DIMS, IMG_DIMS))

        img_count = 0
        for i in range(0, num_images):
            if i + 1 not in bad_frames:
                new_cube[img_count, :, :] = hdu.data[i, :, :]
                img_count += 1

        print("Finished removing %d frames." %num_bad_frames)
        save_cube(new_cube, epoch_name, TRIMMED_SUFFIX)

    except NameError:
        print("Could not load megacube file %s" %cube_path)

if __name__ == "__main__":
    commands = {
        "create"    : create_cube,
        "remove"    : remove_frames,
        "align"     : align_cube
    }

    # sys.argv[1] is command
    # sys.argv[2] is epoch
    command = ""
    path = ""

    try:
        command = sys.argv[1]
        path = sys.argv[2]
    except IndexError:
        print("Error: Invalid number of commands entered")
        print_usage()

    command = command.replace("-", "")
    commands[command](path)