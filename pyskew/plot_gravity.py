#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
# Modified 3/27/2019 by Kevin G.
#############################################

# ***** Import Statements *****
#--------------------------------------------
import os,glob
import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import rasterio
from rasterio.enums import Resampling
import pyproj
import pyskew.utilities as utl
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import re
from functools import cmp_to_key

def get_sandwell(window,down_sample_factor,sandwell_files_path="../raw_data/gravity/Sandwell/*.tiff",resample_method=Resampling.average):
    sandwell_files = glob.glob(sandwell_files_path)

    def cmp_grav_files(x,y): #sorting function utilizing file names to order concatonation of sandwell data
        xfile = os.path.basename(x)
        yfile = os.path.basename(y)
        xlat,xlon = list(map(float,re.findall(r'\d+',xfile)))
        ylat,ylon = list(map(float,re.findall(r'\d+',yfile)))

        if "S" in xfile: xlat = -xlat
        if "W" in xfile: xlon = 360-xlon
        if "S" in yfile: ylat = -ylat
        if "W" in yfile: ylon = 360-ylon

        if xlat-ylat==0: return (utl.convert_to_0_360(ylon)-utl.convert_to_0_360(window[1]))%360 - (utl.convert_to_0_360(xlon)-utl.convert_to_0_360(window[1]))%360
#            if (ylon-window[1])==0: return -1
#            elif (xlon-window[1])==0: return 1
#            else: return (ylon-window[1])%360 - (xlon-window[1])%360
        else: return ylat-xlat

    idx,all_gravs,prev_file,all_lats,all_lons = 0,[],None,np.array([]),np.array([])
    for filepath in sorted(sandwell_files,key=cmp_to_key(cmp_grav_files)):
        with rasterio.open(filepath) as dataset:
            lats = np.arange(dataset.bounds[3],dataset.bounds[1],-dataset.res[1]*down_sample_factor)
            lons = np.arange(dataset.bounds[0],dataset.bounds[2],dataset.res[0]*down_sample_factor)
            if utl.convert_to_0_360(round(window[0],3))>utl.convert_to_0_360(round(dataset.bounds[0],3)): lon_lb = int((utl.convert_to_0_360(window[0])-utl.convert_to_0_360(dataset.bounds[0]))/(dataset.res[0]*down_sample_factor) + .5)
            else: lon_lb = 0
            if utl.convert_to_0_360(round(window[1],3))<utl.convert_to_0_360(round(dataset.bounds[2],3)): lon_ub = int((utl.convert_to_0_360(window[1])-utl.convert_to_0_360(dataset.bounds[2]))/(dataset.res[0]*down_sample_factor) + .5)
            else: lon_ub = -1
            if round(window[2],3)>round(dataset.bounds[1],3): lat_lb = int((window[2]-dataset.bounds[1])/(dataset.res[1]*down_sample_factor) + .5)
            else: lat_lb = 0
            if round(window[3],3)<round(dataset.bounds[3],3): lat_ub = int((window[3]-dataset.bounds[3])/(dataset.res[1]*down_sample_factor) + .5)
            else: lat_ub = -1
            lats = lats[-(lat_ub):-lat_lb-1]
            lons = lons[abs(lon_lb):lon_ub-1]
            if len(lons)==0 or len(lats)==0: continue
            print("Loading: %s"%str(filepath))
            grav = dataset.read(1,
    #            window=gwindow,
                out_shape=(int(dataset.height / down_sample_factor + .5), int(dataset.width / down_sample_factor + .5)),
    #            resampling=Resampling.gauss
    #            resampling=Resampling.nearest
                resampling=resample_method
            )
            grav = grav[-(lat_ub):-lat_lb-1,abs(lon_lb):lon_ub-1]

        if len(lats)>grav.shape[0]: lats = lats[:grav.shape[0]]
        if len(lons)>grav.shape[1]: lons = lons[:grav.shape[1]]

        if prev_file!=None:
            xfile = os.path.basename(prev_file)
            yfile = os.path.basename(filepath)
            xlat,xlon = list(map(float,re.findall(r'\d+',xfile)))
            ylat,ylon = list(map(float,re.findall(r'\d+',yfile)))

            if "S" in xfile: xlat = -xlat
            if "W" in xfile: xlon = 360-xlon
            if "S" in yfile: ylat = -ylat
            if "W" in yfile: ylon = 360-ylon

            if abs(round(xlat-ylat))!=0:
                #create new entry
                idx +=1
                all_lats = np.hstack([all_lats,lats])
                all_gravs.append(grav)
            elif abs(round(xlon-ylon))==60:
                #concatonate left
                if idx==0: all_lons = np.hstack([lons,all_lons])
                all_gravs[idx] = np.hstack([grav,all_gravs[idx]])
            else: raise RuntimeError("Couldn't coordinate how to concatonate %s to gravity array"%str(filepath))
        else:
            all_lats = np.hstack([all_lats,lats])
            all_lons = np.hstack([lons,all_lons])
            all_gravs.append(grav)

        prev_file = filepath

    try:
        all_grav = all_gravs[0]
        for next_grav in all_gravs[1:]:
            all_grav = np.vstack([all_grav,next_grav])
    except (UnboundLocalError,IndexError) as err: return np.array([]),np.array([]),np.array([])

    return all_lons,all_lats,np.array(all_grav)

def plot_deskew_on_grav(*args,**kwargs):
    pass

def plot_deskew_crossings_only(*args,**kwargs):
    pass

def plot_grav(*args,**kwargs):
    pass

def plot_wiggle(*args,**kwargs):
    pass

def plot_chrons_to_map(*args,**kwargs): #maybe
    pass
