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
-cc : float, optional
           sets the convergence criteria for the maximum likelihood pole finder
-ep : Even list of floats, optional
           changes opperation to predict strikes solely from the input euler pole rather than from the data will
           report the degree to which the data predicted strikes agree with the input euler pole as well. Multipule
            Euler poles can be passed and estimates based on all will be reported, however, only the last euler
            pole will be saved in output deskew file

Raises
----------
RuntimeError
ValueError
"""
import os,sys
import numpy as np
import pyrot.max as pymax
from pyrot.rot import cov_to_ellipse,ellipse_to_cov,latlon2cart,cart2latlon
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
from functools import reduce

def calc_strikes_and_add_err(dsk_path,mlat=90,mlon=0,ma=1,mb=1,mphi=0,geoid=Geodesic(6371,0.0),outfile=None,filter_by_quality=False,visualize=False,visual_padding=3.,down_sample_factor=5.,sandwell_files_path="../raw_data/gravity/Sandwell",convergence_level=0.01,euler_pole=None):
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
    convergence_level : float, optional
            the convergence criteria to pass to the function finding the maximum likelihood pole
    euler_pole : iterable, optional
            changes opperation to predict strikes solely from the input euler pole rather than from the data will
            report the degree to which the data predicted strikes agree with the input euler pole as well. Multipule
            Euler poles can be passed and estimates based on all will be reported, however, only the last euler
            pole will be saved in output deskew file

    Returns
    ----------


    Raises
    ----------
    RuntimeError
    ValueError
    """
    (mx,my,mz),mcov = latlon2cart(mlat,mlon,ellipse_to_cov(mlat,mlon,ma,mb,mphi))

    dsk_df = utl.open_deskew_file(dsk_path)
    dsk_df.sort_values("inter_lat",inplace=True,ascending=False)
    if filter_by_quality:
        bad_dsk_data = dsk_df[dsk_df["quality"]!="g"]
        dsk_df = dsk_df[dsk_df["quality"]=="g"]
    tcov,strike_diffs = np.zeros([3,3]),[]
    szs_to_calc = dsk_df["sz_name"].drop_duplicates()#.drop(24) #removes MahiMahi

    if isinstance(euler_pole,type(None)) or len(euler_pole)==0: euler_poles = [None]
    elif len(euler_pole)==2 and (isinstance(euler_pole[0],float) or isinstance(euler_pole[0],int)): euler_poles = [euler_pole]
    elif len(euler_pole)>0 and (isinstance(euler_pole[0],list) or isinstance(euler_pole[0],tuple)): euler_poles = euler_pole
    else: raise ValueError("Euler pole must be None or either a list of euler poles which are length=2 or a single euler pole with lat and lon entries. (i.e. [90,0] or [[90,0],[0,0]])")

    for sz in szs_to_calc:
        sz_df = dsk_df[dsk_df["sz_name"]==sz]
        print(sz,":",len(sz_df.index))

        if visualize:
            window = [utl.convert_to_0_360(sz_df["inter_lon"].min()-visual_padding),utl.convert_to_0_360(sz_df["inter_lon"].max()+visual_padding),sz_df["inter_lat"].min()-visual_padding,sz_df["inter_lat"].max()+visual_padding]
            fig = plt.figure(dpi=100)
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

        num_sites = (sz_df["track_type"]=="aero").sum()/2 + (sz_df["track_type"]=="ship").sum()

        if num_sites>2: #overdetermined case
            data = {"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]}
            for i,row in sz_df.iterrows():
                if row["track_type"]=="aero":
                    if "Ed" in row["comp_name"]: continue
                    elif "Vd" in row["comp_name"]:
                        other_comp = sz_df[sz_df["comp_name"]==row["comp_name"].replace("Vd","Ed")].iloc[0]
                        row["inter_lat"] = (row["inter_lat"]+other_comp["inter_lat"])/2
                        row["inter_lon"] = (row["inter_lon"]+other_comp["inter_lon"])/2
                    else: raise RuntimeError("You really shouldn't have gotten here, you have aeromag that is unrecognized")
                if visualize:
                    if row["quality"]!="g": marker = "X"
                    else:
                        if row["track_type"]=="ship": marker = "o"
                        else: marker = "s"
                    ax.scatter(row["inter_lon"],row["inter_lat"],facecolors=(row["r"],row["g"],row["b"]),edgecolors="k",transform=ccrs.PlateCarree(),marker=marker, zorder=100)
                data["ccl"].append([row["comp_name"],[90.0,0.10,row["inter_lat"],row["inter_lon"]]])
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data,convergence_level=convergence_level)
            for i in range(len(data["ccl"])):
                data["ccl"][i][1][1] *= np.sqrt(chisq)
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data,convergence_level=convergence_level)
            print("\t",(plat,plon,maj_se,min_se,phi),chisq,dof)
            (_,_,_),scov = latlon2cart(plat,plon,ellipse_to_cov(plat,plon,maj_se,min_se,phi))
            tcov += scov
            for ep_idx,euler_pole in enumerate(euler_poles):
                if not isinstance(euler_pole,type(None)):
                    print("--------------------------------------------------------------------------------")
                    print("Euler Pole: %.1f, %.1f"%(euler_pole[0],euler_pole[1]))
                estrikes,dists = [],[]
                for i,row in sz_df.iterrows():
                    if not isinstance(euler_pole,type(None)):
                        geodict = geoid.Inverse(*euler_pole,row["inter_lat"],row["inter_lon"])
                        pgeodict = geoid.Inverse(plat,plon,row["inter_lat"],row["inter_lon"])
                        strike = geodict["azi2"]
                        pstrike = pgeodict["azi2"]+90
                        if pstrike < 0: pstrike += 180
                        strike_diff = abs(strike-pstrike)
                        if strike_diff>90: strike_diff = abs(180-strike_diff)
                        if len(strike_diffs)<ep_idx+1: strike_diffs.append([])
                        strike_diffs[ep_idx].append(strike_diff); estrikes.append(geodict["azi1"]+180)
                    else:
                        pgeodict = geoid.Inverse(plat,plon,row["inter_lat"],row["inter_lon"])
                        strike = pgeodict["azi2"]+90
                    dists.append(pgeodict["a12"])
                    if strike < 0: strike += 180
                    dsk_df.at[i,"strike"] = strike
                    if not isinstance(euler_pole,type(None)): print("\t\t",row["comp_name"],"\n","\t\t\tEuler Pole Strike: ", strike,"\n\t\t\tPredicted Strike: ",pstrike)
                    else: print("\t\t",row["comp_name"], strike)

                if visualize:
                    pdis = np.mean(dists)
                    print("Average Distance to GC Pole: ", pdis)
                    ax = psk.plot_small_circle(plon,plat,pdis,color = "k", m=ax, geoid=Geodesic(6371.,0.), transform=ccrs.PlateCarree(), alpha=.7, linewidth=5, zorder=1)
                    if not isinstance(euler_pole,type(None)):
                        estrike = np.mean(estrikes)
                        print("Average Azimuth of Sites Relative to EP: ", estrike)
                        ep_color = plt.rcParams['axes.prop_cycle'].by_key()['color'][(ep_idx%9)+1]
                        ax = psk.plot_great_circle(euler_pole[1],euler_pole[0],estrike, m=ax, color=ep_color, geoid=Geodesic(6371.,0.), transform=ccrs.PlateCarree(), alpha=.6, linewidth=3, zorder=2)
            if visualize:
                all_lons,all_lats,all_grav = pg.get_sandwell(window,down_sample_factor,resample_method=Resampling.average,sandwell_files_path=os.path.join(sandwell_files_path,"*.tiff"))
                print("Plotting Gravity")
                start_time = time()
                print("Grid Sizes: ",all_lons.shape,all_lats.shape,all_grav.shape)
                fcm = ax.contourf(all_lons, all_lats, all_grav, 60, cmap="Blues_r", alpha=.75, transform=ccrs.PlateCarree(), zorder=0, vmin=0, vmax=255)
                print("Runtime: ",time()-start_time)
                ax.set_extent(window,ccrs.PlateCarree())
                vis_outpath = os.path.join(os.path.dirname(dsk_path),"strike_fit_%s"%sz)
                print("Saving: %s"%vis_outpath)
                fig.savefig(vis_outpath)

        elif num_sites==2: #equal determined case
            strike = geoid.Inverse(sz_df.iloc[0]["inter_lat"],sz_df.iloc[0]["inter_lon"],sz_df.iloc[1]["inter_lat"],sz_df.iloc[1]["inter_lon"])["azi1"]
            if strike < 0: strike += 180
            for i,row in sz_df.iterrows():
                dsk_df.at[i,"strike"] = strike
                print("\t",row["comp_name"], strike)
        else: #under determined case; just ignore
            pass

    if filter_by_quality:
        dsk_df = dsk_df.append(bad_dsk_data)
        dsk_df.sort_values("inter_lat",inplace=True,ascending=False)

    print("--------------------------------------")
    (mlat,mlon),totcov = cart2latlon(mx,my,mz,mcov+tcov)
    full_unc = cov_to_ellipse(mlat,mlon,totcov)
    print("Full Uncertainty: ",full_unc)
    if not isinstance(euler_pole,type(None)):
        if visualize:
            for ep_idx in range(len(strike_diffs)):
                print("For EP %d -> Mean, Median, Min, Max Strike Differences: "%ep_idx,sum(strike_diffs[ep_idx])/len(strike_diffs[ep_idx]),np.median(strike_diffs[ep_idx]),min(strike_diffs[ep_idx]),max(strike_diffs[ep_idx]))
                fig = plt.figure(dpi=100)
                ax = fig.add_subplot(111)
                ax.hist(strike_diffs[ep_idx],bins=np.arange(0.,4.2,0.2))
                ax.axvline(sum(strike_diffs[ep_idx])/len(strike_diffs[ep_idx]),color="tab:red",linestyle="--")
                ax.axvline(np.median(strike_diffs[ep_idx]),color="tab:orange")
                vis_outpath = os.path.join(os.path.dirname(dsk_path),"strike_fit_epstats_%d.png"%ep_idx)
                print("Saving: %s"%vis_outpath)
                fig.savefig(vis_outpath)
            all_strike_diffs = reduce(lambda x,y=[]: x+y, strike_diffs)
            print("For All EP -> Mean, Median, Min, Max Strike Differences: ",sum(all_strike_diffs)/len(all_strike_diffs),np.median(all_strike_diffs),min(all_strike_diffs),max(all_strike_diffs))
            fig = plt.figure(dpi=100)
            ax = fig.add_subplot(111)
            ax.hist(all_strike_diffs,bins=np.arange(0.,4.2,0.2))
            ax.axvline(sum(all_strike_diffs)/len(all_strike_diffs),color="tab:red",linestyle="--")
            ax.axvline(np.median(all_strike_diffs),color="tab:orange")
            vis_outpath = os.path.join(os.path.dirname(dsk_path),"strike_fit_all_epstats.png")
            print("Saving: %s"%vis_outpath)
            fig.savefig(vis_outpath)

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
    if "-cc" in sys.argv:
        kwargs["convergence_level"] = float(sys.argv[sys.argv.index("-cc")+1])
    if "-ep" in sys.argv:
        kwargs["euler_pole"] = []
        for e in sys.argv[sys.argv.index("-ep")+1:]:
            if "-" in e: break
            else: kwargs["euler_pole"].append(float(e))
        if len(kwargs["euler_pole"])<2: kwargs["euler_pole"]=None
        else: kwargs["euler_pole"] = np.array(kwargs["euler_pole"]).reshape([int(len(kwargs["euler_pole"])/2),2]).tolist()
#        kwargs["euler_pole"] = [float(sys.argv[sys.argv.index("-ep")+1]),float(sys.argv[sys.argv.index("-ep")+2])]

    calc_strikes_and_add_err(sys.argv[1],**kwargs)


