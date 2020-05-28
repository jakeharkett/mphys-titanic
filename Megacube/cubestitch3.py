import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

PSF_SOURCE = "psf.fits"
PSF_XY = [95, 111]
IMG_DIMS = 1023
SUFFIX_BADREMOVED = "_badremoved.fits"
SUFFIX_CUBE = ".fits"
CUBE_DATA_OUTPUT_SUFFIX = "stitch_data.csv"

def print_usage():
    print("USAGE:")
    print(" -create FILES_LIST_PATH  Create new cube from file FILES_LIST_PATH")

def resize_img(img):
    return img[:, :IMG_DIMS, :IMG_DIMS]

# Returns a list of paths to input files in .txt file at "path"
def load_target_paths(path):
    target_paths = []

    try:
        target_paths = np.loadtxt(path, dtype=str)
    except:
        raise NameError("Could not find file at %s" %path)

    print("Successfully loaded targets: ")

    for target in target_paths:
        print(target)

    return target_paths

# Counts the number of images in each cube specified in .txt file at "path"
def count_images(path, targets=None):
    target_paths = None

    if targets is None:
        target_paths = load_target_paths(path)
    elif targets is not None:
        target_paths = targets

    img_count = 0

    for target_path in target_paths:
        try:
            hdulist = fits.open(target_path)
            hdu = hdulist[0]
            img_count += hdu.data.shape[0]
        except:
            raise NameError("Path %s not found" %target_path)
    return img_count

def save_cube(cube, cube_name, prefix=None):
    save_name = cube_name + prefix
    hdu_new = fits.PrimaryHDU(cube)
    hdu_new.writeto(save_name, overwrite=True)
    print("Finished writing to %s" %save_name)

def create_cube(target_paths, cube_name):
    cube_data_output_path = cube_name + CUBE_DATA_OUTPUT_SUFFIX
    cube_data_output = None
    try:
        cube_data_output = open(cube_data_output_path, "w")
    except:
        raise NameError("Could not find file %s" %cube_data_output_path)
    
    print("Creating cube %s" %cube_name)
    
    targets = load_target_paths(target_paths)
    num_images = count_images(target_paths, targets)

    megacube = np.zeros((num_images, IMG_DIMS, IMG_DIMS))
    print("Created blank cube of size %d images" %num_images)

    img_count = 0

    print("Beginning stitching of cubes...")

    for index, path in enumerate(targets):
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
        cube_data_output.write(targets[index] + ", " + str(num_images) + ", " + str(img_count) + "\n")

    print("Finished stitching")
    save_cube(megacube, cube_name, SUFFIX_CUBE)
    cube_data_output.close()

def remove_frames(bad_frames_filepath, cube_path):
    bad_frames = []

    try:
        bad_frames = np.loadtxt(bad_frames_filepath, dtype=int)
    except NameError:
        print("Could not load file %s" %bad_frames_filepath)
    except TypeError:
        print("No frames found in badframes file: %s" %bad_frames_filepath)

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
        save_cube(new_cube, cube_path, SUFFIX_BADREMOVED)

    except NameError:
        print("Could not load megacube file %s" %cube_path)


if __name__ == "__main__":
    commands = {
        "create"    : create_cube,
        "remove"    : remove_frames
    }

    # sys.argv[1] is command
    # sys.argv[2] is epoch
    command = ""
    path = ""
    cube_path = ""

    try:
        command = sys.argv[1]
        path = sys.argv[2]
        cube_path = sys.argv[3]
    except IndexError:
        print("Error: Invalid number of commands entered")
        print_usage()

    command = command.replace("-", "")
    commands[command](path, cube_path)
