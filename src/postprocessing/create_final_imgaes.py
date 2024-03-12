import json
import math
import os
import sys
import astropy.cosmology
import astropy.io.fits
import astropy.wcs
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpdaf
import numpy as np
import ref_index
import struct
import project_path_config

mpl.use("TkAgg")


# maybe use PyPDF4
def create_pdf():
	pass


def open_sorted_catalog():
	with open(sys.argv[3], "r") as catalog:
		sorted_catalog = catalog.read()

	return sorted_catalog


def create_final_images():
	sorted_catalog = open_sorted_catalog()

	for line in sorted_catalog:
		pass
