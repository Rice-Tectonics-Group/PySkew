#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
# Modified 3/27/2019 by Kevin G.
#############################################

# ***** Import Statements *****
#--------------------------------------------
#import interpies

import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import rasterio
import rasterio.plot
import pyproj
import pyskew.utilities as util
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def plot_deskew_on_grav(*args,**kwargs):
    pass

def plot_deskew_crossings_only(*args,**kwargs):
    pass

def plot_grav(*args,**kwargs):
    pass

def plot_wiggle(*args,**kwargs,add_crossing=False):
    pass

def plot_chrons_to_map(*args,**kwargs): #maybe
    pass
