import sys
import os
import numpy as np
from astropy.io import fits

if __name__ == "__main__":
    header_files_path = "header_input.txt"

    header_files = np.loadtxt(header_files_path, dtype=str)
    output_file = open("headers_output.csv", "w")

    for header in header_files:
        fileFound = False
        for filename in os.listdir(header):
            if filename.endswith(".fits") and filename.startswith("n0") and not fileFound:
                fileFound = True
                hdulist = fits.open(header + filename, ignore_missing_end=True)
                hdr = hdulist[0].header
                exp_start = hdr["EXPSTART"]
                exp_end = hdr["EXPSTOP"]
                output_file.write(header + "," + exp_start + "," + exp_end + "\n")


    output_file.close()