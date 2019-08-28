#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
#############################################

# ***** Import Statements *****
#--------------------------------------------
#import interpies

import os,sys,glob
import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import utilities as utl
import skewness as sk
from netCDF4 import Dataset as netcdf_dataset
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from time import time
import re
from functools import cmp_to_key

#deskew1 = pd.read_csv('/home/dtw6/Code/PySkew/24r/data/C24r.deskew', sep='\t')
if os.path.isfile(sys.argv[1]): dsk_path = sys.argv[1]
else: raise IOError("File not found: %s"%str(sys.argv[1]))
if "-w" in sys.argv:
    idx = sys.argv.index("-w")
    window = [float(sys.argv[idx+1]),float(sys.argv[idx+2]),float(sys.argv[idx+3]),float(sys.argv[idx+4])] #lon_min,lon_max,lat_min,lat_max
else: window = [160.,240.,-40.,40.]
if "-d" in sys.argv: down_sample_factor = float(sys.argv[sys.argv.index("-d")+1])
else: down_sample_factor = 10
deskew = utl.open_deskew_file(dsk_path)
center_lon = 180
resolution = "10m" #options: 10m, 50m, 110m
landcolor = ""
ax_pos = 111
sandwell_files = glob.glob("../raw_data/gravity/Sandwell/*.tiff")
#sandwell_files = ['../raw_data/gravity/Sandwell/gravN0E180.tiff','../raw_data/gravity/Sandwell/gravN40E180.tiff','../raw_data/gravity/Sandwell/gravN0W120.tiff','../raw_data/gravity/Sandwell/gravN40W120.tiff']
#sandwell_files = ['../raw_data/gravity/Sandwell/gravN40E180.tiff']

# Create map elements
fig = plt.figure(figsize=(16,9),dpi=100)
proj = ccrs.Mercator(central_longitude=center_lon)
ax = fig.add_subplot(ax_pos,projection=proj)
ax.set_xticks(np.arange(0, 370, 10.), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(-80, 90, 10.), crs=ccrs.PlateCarree())
ax.tick_params(grid_linewidth=.5,grid_linestyle=":",color="k",labelsize=8)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
land = cfeature.NaturalEarthFeature('physical', 'land', resolution, edgecolor="black", facecolor=landcolor, linewidth=2)
ax.add_feature(land)
ax.set_extent(window, ccrs.PlateCarree())

def cmp_grav_files(x,y):
    xfile = os.path.basename(x)
    yfile = os.path.basename(y)
    xlat,xlon = list(map(float,re.findall(r'\d+',xfile)))
    ylat,ylon = list(map(float,re.findall(r'\d+',yfile)))

    if "S" in xfile: xlat = -xlat
    if "W" in xfile: xlon = 360-xlon
    if "S" in yfile: ylat = -ylat
    if "W" in yfile: ylon = 360-ylon

    if xlat-ylat==0: return (ylon-window[1])%360 - (xlon-window[1])%360
    else: return ylat-xlat

#for e in sorted(sandwell_files,key=cmp_to_key(cmp_grav_files)):
#    print(e)
#    efile = os.path.basename(e)
#    elat,elon = list(map(float,re.findall(r'\d+',efile)))
#    if "S" in efile: elat = -elat
#    if "W" in efile: elon = 360-elon
#    print(elat,elon)
#    print((elon-window[1])%360)

#sys.exit()

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
idx,all_gravs,prev_file,all_lats,all_lons = 0,[],None,np.array([]),np.array([])
for filepath in sorted(sandwell_files,key=cmp_to_key(cmp_grav_files)):
    print("---------------------------------------------------------------------------------")
    print("Loading: %s"%str(filepath))
    start_time = time()
    with rasterio.open(filepath) as dataset:
        lats = np.arange(dataset.bounds[3],dataset.bounds[1],-dataset.res[1]*down_sample_factor)
        lons = np.arange(dataset.bounds[0],dataset.bounds[2],dataset.res[0]*down_sample_factor)
        print(min(lons),max(lons),min(lats),max(lats),dataset.bounds)
        if window[0]>dataset.bounds[0]: lon_lb = int((window[0]-dataset.bounds[0])/(dataset.res[0]*down_sample_factor) + .5)
        else: lon_lb = 0
        if window[1]<dataset.bounds[2]: lon_ub = int((window[1]-dataset.bounds[2])/(dataset.res[0]*down_sample_factor) + .5)
        else: lon_ub = -1
        if window[2]>dataset.bounds[1]: lat_lb = int((window[2]-dataset.bounds[1])/(dataset.res[1]*down_sample_factor) + .5)
        else: lat_lb = 0
        if window[3]<dataset.bounds[3]: lat_ub = int((window[3]-dataset.bounds[3])/(dataset.res[1]*down_sample_factor) + .5)
        else: lat_ub = -1
        lats = lats[-(lat_ub):-lat_lb-1]
        lons = lons[lon_lb:lon_ub]
        print((int(dataset.height / down_sample_factor), int(dataset.width / down_sample_factor)),dataset.res[1]*down_sample_factor,dataset.res[0]*down_sample_factor)
        print(lon_lb,lon_ub,-(lat_ub),-lat_lb-1)
        print(lons.shape,lats.shape)
        if len(lons)==0 or len(lats)==0: print("No data in window for file %s, skipping"%str(filepath)); continue
        print(min(lons),max(lons),min(lats),max(lats),dataset.bounds)
#        print(lons,lats)
#        gwindow = Window.from_slices(lat_slice,lon_slice, height = int(dataset.height / down_sample_factor), width = int(dataset.width / down_sample_factor))
        grav = dataset.read(1,
#            window=gwindow,
            out_shape=(int(dataset.height / down_sample_factor + .5), int(dataset.width / down_sample_factor + .5)),
#            resampling=Resampling.gauss
#            resampling=Resampling.nearest
            resampling=Resampling.average
        )
        grav = grav[-(lat_ub):-lat_lb-1,lon_lb:lon_ub]

    if len(lats)>grav.shape[0]: lats = lats[:grav.shape[0]]
    if len(lons)>grav.shape[1]: lons = lons[:grav.shape[1]]

    if prev_file!=None:
        xfile = os.path.basename(prev_file)
        yfile = os.path.basename(filepath)
        xlat,xlon = list(map(float,re.findall(r'\d+',xfile)))
        ylat,ylon = list(map(float,re.findall(r'\d+',yfile)))

        if "S" in xfile: xlat = -xlat
        if "W" in xfile: xlon = 180+xlon
        if "S" in yfile: ylat = -ylat
        if "W" in yfile: ylon = 180+ylons

        print(xlat,xlon,ylat,ylon)
        print(round(xlat-ylat),round(xlon-ylon))

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

all_grav = all_gravs[0]
for next_grav in all_gravs[1:]:
    all_grav = np.vstack([all_grav,next_grav])

print(time()-start_time)
print("Plotting")
print(all_lons.shape,all_lats.shape,all_grav.shape)
start_time = time()
fcm = ax.contourf(all_lons, all_lats, all_grav, 60, cmap="viridis", alpha=.75, transform=ccrs.PlateCarree(), zorder=0)
print(time()-start_time)
#cbar = plt.colorbar(fcm)#,ticks=np.arange(0.0,vmax+step,step))

#Older method
#filepath = '../raw_data/gravity/Sandwell/gravN40W120.tiff'
#dataset = rasterio.open(filepath)
#img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
#band1 = dataset.read(1)
#ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)



#for j,(i,row) in enumerate(deskew.iterrows()):
#    # Read in deskewed profile
#    # This is hard-coded now. It will be updated to take a list of profiles in the future
#    dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))

#    # Define the ellipsoid
#    geod = Geodesic.WGS84

#    # Define the angle along which to project
#    perp = row["strike"]-180
#    lon = dskd["lon"]
#    lat = dskd["lat"]
#    mag = sk.phase_shift_data(dskd["mag"].tolist(),row["phase_shift"])

#    # Find distance to project
#    if row["track_type"] == 'ship':
#        pcol = '#000000'
#        scle = 0.2*1e3
#    if row["track_type"] == 'aero':
#        if 'Ed' in row["comp_name"]:
#            pcol = 'purple'
#        else:
#            pcol = 'darkorchid'
#        scle = 0.5*1e3

#    # Project amplitude onto map
#    mlats,mlons = [],[]
#    for i in range(len(mag)):
#        gdsc = geod.Direct(lat[i], lon[i], perp, mag[i]*scle)
#        mlons.append(gdsc['lon2'])
#        mlats.append(gdsc['lat2'])

#    # Plot map elements
#    ax.plot(utl.convert_to_0_360(lon), lat, '--', linewidth=1.0, transform=proj, color=pcol, zorder=990)
#    ax.plot(utl.convert_to_0_360(mlons), mlats, '-', linewidth=1.0, transform=proj, color=pcol, zorder=1000)
#    ax.fill_between(utl.convert_to_0_360(np.array(mlons)[mag>0]), np.array(mlats)[mag>0], lat[mag>0], transform=proj, alpha=0.5, color=pcol)#



#for name,extent in [["All",[240, 200, -15, 35]],["Molokai-Clarion",[225, 220, 16, 23]]]:
#    ax.set_extent(extent, transform=proj)
#    plt.gcf().savefig('/home/dtw6/Documents/bathymetry/' + name + '_bath.png')

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

#ax.add_feature(land, zorder=100)

running,inp = True,"r"
#plt.ion()
plt.show()
#while running:
#    if inp == "r":
#        ax.collections = []
#        deskew = utl.open_deskew_file(dsk_path)
#        for sz_name in deskew["sz_name"].drop_duplicates():
#            ax.scatter(utl.convert_to_0_360(deskew[deskew["sz_name"]==sz_name]["inter_lon"]), deskew[deskew["sz_name"]==sz_name]["inter_lat"], transform=ccrs.PlateCarree(), zorder=10000, s=20, color="k")

#    elif inp == "rt":
#        ax.lines = []
#        deskew = utl.open_deskew_file(dsk_path)
#        for j,(i,row) in enumerate(deskew.iterrows()):
#            # Read in deskewed profile
#            # This is hard-coded now. It will be updated to take a list of profiles in the future
#            dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))

#            # Define the ellipsoid
#            geod = Geodesic.WGS84

#            # Define the angle along which to project
#            perp = row["strike"]-180
#            lon = dskd["lon"]
#            lat = dskd["lat"]
#            mag = sk.phase_shift_data(dskd["mag"].tolist(),row["phase_shift"])

#            # Find distance to project
#            if row["track_type"] == 'ship':
#                pcol = '#000000'
#                scle = 0.2*1e3
#            if row["track_type"] == 'aero':
#                if 'Ed' in row["comp_name"]:
#                    pcol = 'purple'
#                else:
#                    pcol = 'darkorchid'
#                scle = 0.5*1e3

#            # Project amplitude onto map
#            mlats,mlons = [],[]
#            for i in range(len(mag)):
#                gdsc = geod.Direct(lat[i], lon[i], perp, mag[i]*scle)
#                mlons.append(gdsc['lon2'])
#                mlats.append(gdsc['lat2'])

#            # Plot map elements
#            ax.plot(utl.convert_to_0_360(lon), lat, '--', linewidth=1.0, transform=proj, color=pcol, zorder=990)
#            ax.plot(utl.convert_to_0_360(mlons), mlats, '-', linewidth=1.0, transform=proj, color=pcol, zorder=1000)
#            ax.fill_between(utl.convert_to_0_360(np.array(mlons)[mag>0]), np.array(mlats)[mag>0], lat[mag>0], transform=proj, alpha=0.5, color=pcol)#

#    elif inp == "q": running = False; break

#    else: print("Unknown input")

#    plt.draw()
#    plt.pause(.001)
#    inp = input("Plot Command: ")

#plt.close()
