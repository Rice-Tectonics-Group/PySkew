#!/usr/bin/env python3

#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
# Modified 3/27/2019 by Kevin G.
#############################################

# ***** Import Statements *****
#--------------------------------------------
import os,sys
import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyskew.utilities as utl
import pyskew.skewness as sk
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import pyskew.plot_gravity as pg
from time import time
from rasterio.enums import Resampling

#deskew1 = pd.read_csv('/home/dtw6/Code/PySkew/24r/data/C24r.deskew', sep='\t')
if os.path.isfile(sys.argv[1]): dsk_path = sys.argv[1]
else: raise IOError("File not found: %s"%str(sys.argv[1]))
if "-w" in sys.argv:
    idx = sys.argv.index("-w")
    window = [utl.convert_to_0_360(sys.argv[idx+1]),utl.convert_to_0_360(sys.argv[idx+2]),float(sys.argv[idx+3]),float(sys.argv[idx+4])] #lon_min,lon_max,lat_min,lat_max
else: window = [160.,240.,-40.,40.]
if "-d" in sys.argv: down_sample_factor = float(sys.argv[sys.argv.index("-d")+1])
else: down_sample_factor = 10
if "-swf" in sys.argv: sandwell_files_path = sys.argv[sys.argv.index("-swf")+1]
else: sandwell_files_path = "../raw_data/gravity/Sandwell"
if "-rb" in sys.argv: remove_bad_data = True
else: remove_bad_data = False
center_lon = 180
resolution = "10m" #options: 10m, 50m, 110m
landcolor = ""
ax_pos = 111
# Define the ellipsoid
geod = Geodesic(6371000.,0.0) #Radius in meters

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

## Create map elements
#fig = plt.figure(figsize=(16,9),dpi=100)
#proj = ccrs.Mollweide()
#ax = fig.add_subplot(111,projection=proj)
#plt.gca().set_axis_off()
#plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
#plt.margins(0,0)
#land = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='black', facecolor='gray', linewidth=0.5,zorder=1000)

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
#cbar = plt.colorbar(fcm)#,ticks=np.arange(0.0,vmax+step,step))

#Older method
#filepath = '../raw_data/gravity/Sandwell/gravN40W120.tiff'
#dataset = rasterio.open(filepath)
#img_extent = (dataset.bounds[0], dataset.bounds[2], dataset.bounds[1], dataset.bounds[3])
#band1 = dataset.read(1)
#ax.imshow(band1, origin='upper', extent=img_extent, transform=proj, zorder=0, alpha=0.75)

running,inp = True,"r"
plt.ion()
plt.show()
while running:

    if "rl" in inp or inp.split()[0] == "r":
        try:
            for sp in site_points:
                sp.remove()
        except NameError: pass
        deskew = utl.open_deskew_file(dsk_path)
        print("Plotting Sites")
        site_points = []
        for i,sz_name in enumerate(deskew["sz_name"].drop_duplicates()):
            sz_dsk = deskew[deskew["sz_name"]==sz_name]
            bad_data = sz_dsk[sz_dsk["quality"]!="g"]
            bad_color = bad_data[["r","g","b"]].to_numpy()
            good_data = sz_dsk[sz_dsk["quality"]=="g"]
            aero_df = good_data[good_data["track_type"]=="aero"]
            ship_df = good_data[good_data["track_type"]=="ship"]
            aero_color = aero_df[["r","g","b"]].to_numpy()
            ship_color = ship_df[["r","g","b"]].to_numpy()
#            color = (sz_dsk.iloc[0]["r"],sz_dsk.iloc[0]["g"],sz_dsk.iloc[0]["b"])
            site_points.append(ax.scatter(utl.convert_to_0_360(bad_data["inter_lon"]), bad_data["inter_lat"], transform=ccrs.PlateCarree(), zorder=10000, s=20, color=bad_color, edgecolor="k", marker="X"))
            site_points.append(ax.scatter(utl.convert_to_0_360(aero_df["inter_lon"]), aero_df["inter_lat"], transform=ccrs.PlateCarree(), zorder=10000, s=20, color=aero_color, edgecolor="k", marker="s"))
            site_points.append(ax.scatter(utl.convert_to_0_360(ship_df["inter_lon"]), ship_df["inter_lat"], transform=ccrs.PlateCarree(), zorder=10000, s=20, color=ship_color, edgecolor="k", marker="s"))

    if "rt" in inp or inp.split()[0] == "r":
        try:
            for dt in deskew_tracks:
                dt[0].remove()
            for df in deskew_fill:
                df[0].remove()
        except NameError: pass
        except TypeError: pass
        print("Plotting Tracks: ")
        deskew,deskew_tracks,deskew_fill = utl.open_deskew_file(dsk_path),[],[]
        for j,(i,row) in enumerate(deskew.iterrows()):
            if remove_bad_data and row["quality"]!="g": continue
            # Read in deskewed profile
            # This is hard-coded now. It will be updated to take a list of profiles in the future
            dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))
            print("\tPlotting: ", row["comp_name"])

            # Define the angle along which to project
            perp = row["strike"]-180
            lon = dskd["lon"].tolist()
            lat = dskd["lat"].tolist()
            mag = sk.phase_shift_data(dskd["mag"].tolist(),row["phase_shift"])
#            print("\t\t",perp)
#            print("\t\t",np.array(lon).shape,np.array(lat).shape,np.array(mag).shape)

            # Find distance to project
            if row["track_type"] == 'ship':
                pcol = '#000000'
                scle = 0.2*1e3
            elif row["track_type"] == 'aero':
                if 'Ed' in row["comp_name"]:
                    pcol = 'purple'
                else:
                    pcol = 'darkorchid'
                scle = 0.5*1e3
            else: raise ValueError("Apparently neither an airplane or ship hauled this special magnetometer check your track_type")

            if row["quality"]!="g":
                if row["track_type"] == 'ship': pcol = "lightgrey"
                elif row["track_type"] == 'aero': pcol = "plum"
                alpha,fill_alpha = .5,.2
            else:
                alpha,fill_alpha = 1.,.5

            # Project amplitude onto map
            mlats,mlons,fill_lats = [],[],[]
            for i in range(len(mag)):
                gdsc = geod.Direct(lat[i], lon[i], perp, mag[i]*scle)
                mlons.append(gdsc['lon2'])
                mlats.append(gdsc['lat2'])
                if mag[i]>0: fill_lats.append(gdsc['lat2'])
                else: fill_lats.append(gdsc['lat1'])

            # Plot map elements
            geotrack = proj.transform_points(ccrs.PlateCarree(),np.array(utl.convert_to_0_360(lon)),np.array(lat))
            magtrack = proj.transform_points(ccrs.PlateCarree(),np.array(utl.convert_to_0_360(mlons)),np.array(mlats))
            filltrack = proj.transform_points(ccrs.PlateCarree(),np.array(utl.convert_to_0_360(mlons)),np.array(fill_lats))
            deskew_tracks.append(ax.plot(geotrack[:,0], geotrack[:,1], '--', linewidth=1.0, color=pcol, zorder=990, alpha=alpha))
            deskew_tracks.append(ax.plot(magtrack[:,0], magtrack[:,1], '-', linewidth=1.0, color=pcol, zorder=1000, alpha=alpha))
            deskew_fill.append(ax.fill_between(filltrack[:,0], filltrack[:,1], geotrack[:,1], alpha=fill_alpha, color=pcol))

    if "rg" in inp or inp.split()[0] == "r":

        inpl = inp.split()

        if "-w" in inpl:
            idx = inpl.index("-w")
            window = [utl.convert_to_0_360(inpl[idx+1]),utl.convert_to_0_360(inpl[idx+2]),float(inpl[idx+3]),float(inpl[idx+4])] #lon_min,lon_max,lat_min,lat_max
        if "-d" in inpl: down_sample_factor = float(inpl[inpl.index("-d")+1])

        print("Grav Window and Down sample: ",window,down_sample_factor)
        all_lons,all_lats,all_grav = pg.get_sandwell(window,down_sample_factor,resample_method=Resampling.average,sandwell_files_path=os.path.join(sandwell_files_path,"*.tiff"))

        print("Plotting Gravity")
        start_time = time()
        print("Grid Sizes: ",all_lons.shape,all_lats.shape,all_grav.shape)
        fcm = ax.contourf(all_lons, all_lats, all_grav, 51, cmap="Blues_r", alpha=.75, transform=ccrs.PlateCarree(), zorder=0, vmin=0, vmax=255)
        print("Runtime: ",time()-start_time)

        ax.set_extent(window, ccrs.PlateCarree())

    if inp == "q": running = False; break

    ax.set_facecolor('w')
    plt.draw()
    plt.pause(.001)
    inp = input("Plot Command: ")

plt.close()
