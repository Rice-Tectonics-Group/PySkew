#############################################
# A script to read, manipulate, and test
# the data published in Konrad et al. (2018)
#
# Created 22/01/2019 by Daniel W.
#############################################

# ***** Import Statements *****
#--------------------------------------------
#import interpies

import os,sys
import pandas as pd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import pyproj
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import utilities as utl
import skewness as sk
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import pyskew.plot_gravity as pg
from time import time

#deskew1 = pd.read_csv('/home/dtw6/Code/PySkew/24r/data/C24r.deskew', sep='\t')
if os.path.isfile(sys.argv[1]): dsk_path = sys.argv[1]
else: raise IOError("File not found: %s"%str(sys.argv[1]))
if "-w" in sys.argv:
    idx = sys.argv.index("-w")
    window = [float(sys.argv[idx+1]),float(sys.argv[idx+2]),float(sys.argv[idx+3]),float(sys.argv[idx+4])] #lon_min,lon_max,lat_min,lat_max
else: window = [160.,240.,-40.,40.]
if "-d" in sys.argv: down_sample_factor = float(sys.argv[sys.argv.index("-d")+1])
else: down_sample_factor = 10
center_lon = 180
resolution = "10m" #options: 10m, 50m, 110m
landcolor = ""
ax_pos = 111

# Create map elements
fig = plt.figure(figsize=(9,9),dpi=100)
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

    if inp == "rl" or inp == "r":
        try:
            for sp in site_points:
                ax.collection.remove(sp)
        except NameError: pass
        deskew = utl.open_deskew_file(dsk_path)
        print("Plotting Sites")
        for sz_name in deskew["sz_name"].drop_duplicates():
            site_points = ax.scatter(utl.convert_to_0_360(deskew[deskew["sz_name"]==sz_name]["inter_lon"]), deskew[deskew["sz_name"]==sz_name]["inter_lat"], transform=ccrs.PlateCarree(), zorder=10000, s=20, color="k")

    if inp == "rt" or inp == "r":
        try:
            for dt in deskew_tracks:
                ax.lines.remove(dt)
            for df in deskew_fill:
                ax.collection.remove(df)
        except NameError: pass
        print("Plotting Tracks: ")
        deskew,deskew_tracks,deskew_fill = utl.open_deskew_file(dsk_path),[],[]
        for j,(i,row) in enumerate(deskew.iterrows()):
            # Read in deskewed profile
            # This is hard-coded now. It will be updated to take a list of profiles in the future
            dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))
            print("\tPlotting: ", row["comp_name"])

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
            deskew_tracks.append(ax.plot(utl.convert_to_0_360(lon), lat, '--', linewidth=1.0, transform=ccrs.PlateCarree(), color=pcol, zorder=990))
            deskew_tracks.append(ax.plot(utl.convert_to_0_360(mlons), mlats, '-', linewidth=1.0, transform=ccrs.PlateCarree(), color=pcol, zorder=1000))
            deskew_fill.append(ax.fill_between(utl.convert_to_0_360(np.array(mlons)[mag>0]), np.array(mlats)[mag>0], lat[mag>0], transform=ccrs.PlateCarree(), alpha=0.5, color=pcol))

    if inp == "rg" or inp == "r":

        inpl = inp.split()

        if "-w" in inpl:
            idx = inpl.index("-w")
            window = [float(inpl[idx+1]),float(inpl[idx+2]),float(inpl[idx+3]),float(inpl[idx+4])] #lon_min,lon_max,lat_min,lat_max
        if "-d" in inpl: down_sample_factor = float(inpl[inpl.index("-d")+1])

        all_lons,all_lats,all_grav = pg.get_sandwell(window,down_sample_factor)

        print("Plotting Gravity")
        start_time = time()
        print(all_lons.shape,all_lats.shape,all_grav.shape)
        fcm = ax.contourf(all_lons, all_lats, all_grav, 60, cmap="viridis", alpha=.75, transform=ccrs.PlateCarree(), zorder=0)
        print("Runtime: ",time()-start_time)

        ax.set_extent(window, ccrs.PlateCarree())

    if inp == "q": running = False; break

    plt.draw()
    plt.pause(.001)
    inp = input("Plot Command: ")

plt.close()
