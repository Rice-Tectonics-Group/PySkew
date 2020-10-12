"""
Script to quickly recalculate strikes from anomaly crossings identified using SynthMagGUI
and the greater PySkew library. Strikes are calculated by inverting for a maximum likelihood
great circle using a python version of the "max" program originally written by Richard Gordon.
For more information on the algorithm see Gordon & Cox 1980. Requires the seperate PyRot
library to be installed in PYTHONPATH env variable so it can use its implementation of "max".
If a paleomagnetic pole calculated from these anomalies is provided via the -ell flag then
the average strike uncertainty calculated by max will be added to the pole for a final result.
Note: the function this script runs prints excessively as a diagnostic because many things
can be easily misinterpreted about this function, if this is undesireable pipe to /dev/null.

Parameters
----------
dsk_path : str
           Path to deskew file containing the profiles to invert strikes on

Flags
----------
-h : no input
           Prints this help message and quits
-ell : 5 floats, optional
           Paleomagnetic pole to add strike uncertainty to
-o : str, optional
           Outfile to save new deskew file, if omitted is saved as $INFILE_strike_cor.deskew
-fp : no input
           filters out bad data so strikes are fit only to uncommented records
-v : no input
           creates images visualizing the fit to the sites
-fg : str, optional
           allows user to provide the path to the .tiff gravity files to render under fit

Raises
----------
RuntimeError
"""
import os,sys
import numpy as np
import pyrot.max as pymax
from pyrot.rot import cov_to_ellipse,ellipse_to_cov
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import pyskew.plot_skewness as psk
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import pyskew.plot_gravity as pg
from time import time
from rasterio.enums import Resampling

def calc_strikes_and_add_err(dsk_path,mlat=90,mlon=0,ma=1,mb=1,mphi=0,geoid=Geodesic(6371,0.0),outfile=None,filter_by_quality=False,visualize=False,visual_padding=3.,down_sample_factor=5.,sandwell_files_path="../raw_data/gravity/Sandwell"):
    """
    Function that does the heavy lifting calculating the great circles and associated strikes
    for anomaly crossings. Will also add average strike uncertainty to a paleomagnetic pole
    if provided. Function prints quite a lot for user diagnostics, pipe to /dev/null to silence.

    Parameters
    ----------
    dsk_path : str
            Path to deskew file containing the profiles to invert strikes on
    mlat : float, optional
    mlon : float, optional
    ma : float, optional
    mb : float, optional
    mphi : float, optional
    geoid : geographiclib.geodesic.Geodesic, optional
            geodesic to use for projection default is sphere radius 6371 flattening of 0
    outfile : str, optional
            output deskew file with the corrected strikes
    filter_by_quality : bool, optional
            bool that descides if "bad" data is filtered out of strike fit (Default : False)
    visualize : bool, optional
            weather or not you want to render images showing the fit to sites (Default : False)
    visual_padding : float, optional
            how much to pad out plotting window in degrees (Default : 3.)
    down_sample_factor : float, optional
            how much to downsample the gravity data operates as a divisor so 2 means half of the gravity will plot
            (Default : 5.)
    sandwell_files_path : str, optional
            path to the files containing the sandwell gravity grid to render in .tiff format

    Returns
    ----------


    Raises
    ----------
    RuntimeError
    """
    mcov = ellipse_to_cov(mlat,mlon,ma,mb,mphi)

    dsk_df = utl.open_deskew_file(dsk_path)
    dsk_df.sort_values("inter_lat",inplace=True,ascending=False)
    if filter_by_quality:
        bad_dsk_data = dsk_df[dsk_df["quality"]!="g"]
        dsk_df = dsk_df[dsk_df["quality"]=="g"]
    tcov = np.zeros([2,2])
    szs_to_calc = dsk_df["sz_name"].drop_duplicates()#.drop(24) #removes MahiMahi

    for sz in szs_to_calc:
        sz_df = dsk_df[dsk_df["sz_name"]==sz]
        print(sz,":",len(sz_df.index))

        if visualize:
            window = [utl.convert_to_0_360(sz_df["inter_lon"].min()-visual_padding),utl.convert_to_0_360(sz_df["inter_lon"].max()+visual_padding),sz_df["inter_lat"].min()-visual_padding,sz_df["inter_lat"].max()+visual_padding]
            fig = plt.figure(figsize=(16,9),dpi=100)
            proj = ccrs.Mercator(central_longitude=sz_df["inter_lon"].mean())
            ax = fig.add_subplot(111,projection=proj)
            ax.set_xticks(np.arange(0, 370, 10.), crs=ccrs.PlateCarree())
            ax.set_yticks(np.arange(-80, 90, 10.), crs=ccrs.PlateCarree())
            ax.tick_params(grid_linewidth=.5,grid_linestyle=":",color="k",labelsize=8)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            land = cfeature.NaturalEarthFeature('physical', 'land', "50m", edgecolor="black", facecolor="grey", linewidth=2)
            ax.add_feature(land)

        if len(sz_df.index)>2: #overdetermined case
            data = {"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]}
            for i,row in sz_df.iterrows():
                if row["track_type"]=="aero":
                    if "Ed" in row["comp_name"]: continue
                    elif "Vd" in row["comp_name"]:
                        other_comp = sz_df[sz_df["comp_name"]==row["comp_name"].replace("Vd","Ed")].iloc[0]
                        row["inter_lat"] = (row["inter_lat"]+other_comp["inter_lat"])/2
                        row["inter_lon"] = (row["inter_lon"]+other_comp["inter_lon"])/2
                    else: raise RuntimeError("You really shouldn't have gotten here, you have aeromag that is unrecognized")
                if visualize: ax.scatter(row["inter_lon"],row["inter_lat"],facecolors=(row["r"],row["g"],row["b"]),edgecolors="k",transform=ccrs.PlateCarree())
                data["ccl"].append([row["comp_name"],[90.0,0.10,row["inter_lat"],row["inter_lon"]]])
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data)
            for i in range(len(data["ccl"])):
                data["ccl"][i][1][1] *= np.sqrt(chisq)
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data)
            print("\t",(plat,plon,maj_se,min_se,phi),chisq,dof)
            scov = ellipse_to_cov(plat,plon,maj_se,min_se,phi)
            tcov += scov
            dists = []
            for i,row in sz_df.iterrows():
                geodict = geoid.Inverse(plat,plon,row["inter_lat"],row["inter_lon"])
                strike = geodict["azi2"]+90
                dists.append(geodict["a12"])
                if strike < 0: strike += 180
                dsk_df.at[i,"strike"] = strike
                print("\t\t",row["comp_name"], strike)

            if visualize:
                pdis = np.mean(dists)
                print("Average Distance to Pole: ", pdis)
                ax = psk.plot_small_circle(plon,plat,pdis,color = "k",m=ax,geoid=Geodesic(6371.,0.),transform=ccrs.PlateCarree())
                all_lons,all_lats,all_grav = pg.get_sandwell(window,down_sample_factor,resample_method=Resampling.average,sandwell_files_path=os.path.join(sandwell_files_path,"*.tiff"))
                print("Plotting Gravity")
                start_time = time()
                print("Grid Sizes: ",all_lons.shape,all_lats.shape,all_grav.shape)
                fcm = ax.contourf(all_lons, all_lats, all_grav, 60, cmap="Blues_r", alpha=.75, transform=ccrs.PlateCarree(), zorder=0, vmin=0, vmax=255)
                print("Runtime: ",time()-start_time)
                ax.set_extent(window,ccrs.PlateCarree())
                fig.savefig("strike_fit_%s"%sz)

        else: #under or equal determined case
            strike = geoid.Inverse(sz_df.iloc[0]["inter_lat"],sz_df.iloc[0]["inter_lon"],sz_df.iloc[1]["inter_lat"],sz_df.iloc[1]["inter_lon"])["azi1"]
            if strike < 0: strike += 180
            for i,row in sz_df.iterrows():
                dsk_df.at[i,"strike"] = strike
                print("\t",row["comp_name"], strike)

    if filter_by_quality:
        dsk_df = dsk_df.append(bad_dsk_data)
        dsk_df.sort_values("inter_lat",inplace=True,ascending=False)

    print("--------------------------------------")
    full_unc = cov_to_ellipse(mlat,mlon,mcov+tcov/len(szs_to_calc))
    print("Full Uncertainty: ",full_unc)

    if isinstance(outfile,type(None)): outfile = os.path.join(os.path.dirname(dsk_path),"strike_cor_"+os.path.basename(dsk_path))
    print("Writing to %s"%str(outfile))
    utl.write_deskew_file(outfile,dsk_df)

    return full_unc

if __name__=="__main__":

    kwargs = {}
    if "-h" in sys.argv:
        help(__name__)
        sys.exit()
    if "-ell" in sys.argv:
        iell = sys.argv.index("-ell")
        kwargs["mlat"] = float(sys.argv[iell+1])
        kwargs["mlon"] = float(sys.argv[iell+2])
        kwargs["ma"] = float(sys.argv[iell+3])
        kwargs["mb"] = float(sys.argv[iell+4])
        kwargs["mphi"] = float(sys.argv[iell+5])
    if "-o" in sys.argv:
        io = sys.argv.index("-o")
        kwargs["outfile"] = sys.argv[io+1]
    if "-fq" in sys.argv:
        kwargs["filter_by_quality"] = True
    if "-v" in sys.argv:
        kwargs["visualize"] = True
    if "-fg" in sys.argv:
        kwargs["sandwell_files_path"] = sys.argv[sys.argv.index("-fg")+1]

    calc_strikes_and_add_err(sys.argv[1],**kwargs)
