#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
#############################################

# ***** Import Statements *****
#--------------------------------------------
#import interpies

import os
import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import rasterio
import rasterio.plot
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import utilities as utl
import skewness as sk

#deskew1 = pd.read_csv('/home/dtw6/Code/PySkew/24r/data/C24r.deskew', sep='\t')
deskew = utl.open_deskew_file('../27r/data/pole_9_78.deskew')
#deskew = pd.concat([deskew1,deskew2])
#deskew.reset_index(drop=True, inplace=True)
#cross_lat = deskew['inter_lat']
#cross_lon = deskew['inter_lon']
#deskew_files = deskew['comp_name']
#mag_type = deskew['track_type']

# Create map elements
fig = plt.figure(figsize=(16,9),dpi=100)
proj = ccrs.Mollweide()
ax = fig.add_subplot(111,projection=proj)
plt.gca().set_axis_off()
plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
plt.margins(0,0)
land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5,zorder=1000)

## Read Barckhausen data
#ccz_chrons, gcz_chrons = utl.get_barckhausen_2013_chrons(barckhausen_path='../raw_data/chrons/Barckhausen2013/GSFML.Barckhausen++_2013_MGR.picks.gmt')

#for chron,v in ccz_chrons.items():
#    isochron = ccz_chrons[chron]
#    cczdf = pd.DataFrame(isochron)
#    if not cczdf.empty:
#        ax.plot(utl.convert_to_0_360(cczdf[0]),cczdf[1], transform=proj, zorder=9000, color='gray', linewidth=1)
#for chron,v in gcz_chrons.items():
#    isochron = gcz_chrons[chron]
#    gczdf = pd.DataFrame(isochron)
#    if not gczdf.empty:
#        ax.plot(utl.convert_to_0_360(gczdf[0]),gczdf[1], transform=proj, zorder=9000, color='gray', linewidth=1)

# Plot
filepath = '../raw_data/gravity/Sandwell/gravN40E180.tiff'
dataset = rasterio.open(filepath)
img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
band1 = dataset.read(1)
ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)

filepath = '../raw_data/gravity/Sandwell/gravN0E180.tiff'
dataset = rasterio.open(filepath)
img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
band1 = dataset.read(1)
ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)

filepath = '../raw_data/gravity/Sandwell/gravN0W120.tiff'
dataset = rasterio.open(filepath)
img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
band1 = dataset.read(1)
ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)

filepath = '../raw_data/gravity/Sandwell/gravN40W120.tiff'
dataset = rasterio.open(filepath)
img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
band1 = dataset.read(1)
ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)

for sz_name in deskew["sz_name"].drop_duplicates():
    ax.scatter(utl.convert_to_0_360(deskew[deskew["sz_name"]==sz_name]["inter_lon"]), deskew[deskew["sz_name"]==sz_name]["inter_lat"], transform=proj, zorder=10000, s=20)

## Plot Lin's crossings
#lon20r = [222.818, 222.143, 222.498, 222.742]
#lat20r = [15.288, 12.755, 11.113, 10.136]
#ax.scatter(utl.convert_to_0_360(lon20r), lat20r, transform=proj, zorder=10000, color='red', s=20)

## Plot Lily's crossings
#lon26r = [215.873, 215.990, 215.554, 216.402, 216.451, 216.850, 216.815]
#lat26r = [12.896,  12.438,  11.065,  8.424,   8.006,   6.343,   6.189]
#ax.scatter(utl.convert_to_0_360(lon26r), lat26r, transform=proj, zorder=10000, color='blue', s=20)

## Plot Katerina's crossings
#lon25r = [216.9, 217.6, 217.7]
#lat25r = [8.4, 6.2, 6.1]
#ax.scatter(utl.convert_to_0_360(lon25r), lat25r, transform=proj, zorder=10000, color='green', s=20)

ax.add_feature(land, zorder=100)

for j,(i,row) in enumerate(deskew.iterrows()):
    # Read in deskewed profile
    # This is hard-coded now. It will be updated to take a list of profiles in the future
    dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))

    # Define the ellipsoid
    geod = Geodesic.WGS84

    # Define the angle along which to project
    perp = row["strike"]-180
    lon = dskd["lon"]
    lat = dskd["lat"]
    mag = sk.phase_shift_data(dskd["mag"].tolist(),row["phase_shift"])

    # Find distance to project
    if row["track_type"] == 'ship':
        pcol = '#000000'
        scle = 0.2*1e3
    if row["track_type"] == 'aero':
        if 'Ed' in row["comp_name"]:
            pcol = 'purple'
        else:
            pcol = 'darkorchid'
        scle = 0.5*1e3

    # Project amplitude onto map
    mlats,mlons = [],[]
    for i in range(len(mag)):
        gdsc = geod.Direct(lat[i], lon[i], perp, mag[i]*scle)
        mlons.append(gdsc['lon2'])
        mlats.append(gdsc['lat2'])

    # Plot map elements
    ax.plot(utl.convert_to_0_360(lon), lat, '--', linewidth=1.0, transform=proj, color=pcol, zorder=990)
    ax.plot(utl.convert_to_0_360(mlons), mlats, '-', linewidth=1.0, transform=proj, color=pcol, zorder=1000)
    ax.fill_between(utl.convert_to_0_360(np.array(mlons)[mag>0]), np.array(mlats)[mag>0], lat[mag>0], transform=proj, alpha=0.5, color=pcol)#

#for name,extent in [["All",[240, 200, -15, 35]],["Molokai-Clarion",[225, 220, 16, 23]]]:
#    ax.set_extent(extent, transform=proj)
#    plt.gcf().savefig('/home/dtw6/Documents/bathymetry/' + name + '_bath.png')

plt.show()
