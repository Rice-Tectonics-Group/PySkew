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
-max : Max file to use in calculating strike uncertainty corrected paleomagnetic uncertainties
-askw : Flag which causes max to treat anomalous skewness as an adjustable parameter during paleomag uncertainty calculation
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

def calc_strikes_and_add_err(dsk_path,max_file=None,solve_anom_skew=False,geoid=Geodesic(6371,0.0),outfile=None,filter_by_quality=False,visualize=False,visual_padding=3.,down_sample_factor=5.,sandwell_files_path="../raw_data/gravity/Sandwell",convergence_level=0.01,euler_pole=None):
    """
    Function that does the heavy lifting calculating the great circles and associated strikes
    for anomaly crossings. Will also add average strike uncertainty to a paleomagnetic pole
    if provided. Function prints quite a lot for user diagnostics, pipe to /dev/null to silence.

    Parameters
    ----------
    dsk_path : str
            Path to deskew file containing the profiles to invert strikes on
    max_file : str
            Path to a max file which can be read to calculate a paleomagnetic pole and to add the
            strike uncertainty to the paleomagnetic uncertainty
    solve_anom_skew : bool
            determines if anomalous skewness should be an adjustable parameter in the paleomagnetic
            pole calculation done by max
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
    #First deal with determining the paleomagnetic result if we're updating it's uncertainties
    if max_file!=None:
        comment,header,pdata = pymax.read_max_file(".tmp.max")
        if len(pdata["phs"])==0: raise ValueError("Empty Max file read during correction for strike uncertainty")
        if solve_anom_skew: (mlat,mlon,mmag,askw,ma,mb,mphi),chisq,dof = pymax.max_likelihood_pole(pdata, trial_pole=header[:3], out_path="calc_strike_before.maxout", step=header[-1], max_steps=100, comment=comment, solve_anom_skew=solve_anom_skew); print("Anomalous Skewness Estimated")
        else: (mlat,mlon,mmag,ma,mb,mphi),chisq,dof = pymax.max_likelihood_pole(pdata, trial_pole=header[:3], out_path="calc_strike_before.maxout", step=header[-1], max_steps=100, comment=comment, solve_anom_skew=solve_anom_skew)
        (mx,my,mz),mcov = latlon2cart(mlat,mlon,ellipse_to_cov(mlat,mlon,ma,mb,mphi))
        print("Original Pole: ", mlat,mlon,ma,mb,mphi,chisq,dof)
        print("Original Covariance Matrix:\n", mcov)
        #Make empty preturbed data dictionaries
        pdata_splus,pdata_sminus = {"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]},{"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]}

    #Then get the deskew data we need to estimate strikes from and initalize a few variables
    dsk_df = utl.open_deskew_file(dsk_path)
    dsk_df.sort_values("inter_lat",inplace=True,ascending=False)
    if filter_by_quality: #If removing commented data, remove it
        bad_dsk_data = dsk_df[dsk_df["quality"]!="g"]
        dsk_df = dsk_df[dsk_df["quality"]=="g"]
    strike_diffs = []
    szs_to_calc = dsk_df["sz_name"].drop_duplicates()#.drop(24) #removes MahiMahi

    #Deal with Euler poles if any have been provided
    if isinstance(euler_pole,type(None)) or len(euler_pole)==0: euler_poles = [None]
    elif len(euler_pole)==2 and (isinstance(euler_pole[0],float) or isinstance(euler_pole[0],int)): euler_poles = [euler_pole]
    elif len(euler_pole)>0 and (isinstance(euler_pole[0],list) or isinstance(euler_pole[0],tuple)): euler_poles = euler_pole
    else: raise ValueError("Euler pole must be None or either a list of euler poles which are length=2 or a single euler pole with lat and lon entries. (i.e. [90,0] or [[90,0],[0,0]])")

    n = 0 #record the number of szs strikes are estimated from just in case
    for sz in szs_to_calc: #Iterate over spreading zones
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

            #Construct Max file for Strike Estimation
            for i,row in sz_df.iterrows():
                if row["track_type"]=="aero":
                    if "Ed" in row["comp_name"]: continue
                    elif "Vd" in row["comp_name"]:
                        other_comp = sz_df[sz_df["comp_name"]==row["comp_name"].replace("Vd","Ed")].iloc[0]
                        row["inter_lat"] = (row["inter_lat"]+other_comp["inter_lat"])/2
                        row["inter_lon"] = (row["inter_lon"]+other_comp["inter_lon"])/2
                    else: raise RuntimeError("You really shouldn't have gotten here, you have aeromag that can't find its second component")
                #Visualize Estimate against bathymetery and sites
                if visualize:
                    if row["quality"]!="g": marker = "X"
                    else:
                        if row["track_type"]=="ship": marker = "o"
                        else: marker = "s"
                    ax.scatter(row["inter_lon"],row["inter_lat"],facecolors=(row["r"],row["g"],row["b"]),edgecolors="k",transform=ccrs.PlateCarree(),marker=marker, zorder=100)
                data["ccl"].append([row["comp_name"],[90.0,0.10,row["inter_lat"],row["inter_lon"]]])
            #Estimate Strike With Max, correct for over/underdispursion
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data,convergence_level=convergence_level)
            for i in range(len(data["ccl"])):
                data["ccl"][i][1][1] *= np.sqrt(chisq)
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data,convergence_level=convergence_level)
            print("\t",(plat,plon,maj_se,min_se,phi),chisq,dof)

            if max_file!=None: #Create Preturbed Max data for Strike Uncertainty
                for i,row in sz_df.iterrows(): #Find all max data entries associated with this SZ
                    pdata_idx = find_pdata_idx(row["comp_name"], pdata)
                    if pdata_idx==None: continue #This is happening because of aeromag
                    pdata_eplus = [pdata["phs"][pdata_idx][0],pdata["phs"][pdata_idx][1].copy()]
                    pdata_eplus[1][-1] += maj_se #Add strike SE
                    pdata_splus["phs"].append(pdata_eplus)
                    pdata_eminus = [pdata["phs"][pdata_idx][0],pdata["phs"][pdata_idx][1].copy()]
                    pdata_eminus[1][-1] -= maj_se #Subtract strike SE
                    pdata_sminus["phs"].append(pdata_eminus)
#            (_,_,_),scov = latlon2cart(plat,plon,ellipse_to_cov(plat,plon,maj_se,min_se,phi))
#            tcov += scov
            n += 1

            #If Estimating from Euler poles instead
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
                    if strike < 0: strike += 360
                    if strike < 180: strike += 180
                    dsk_df.at[i,"strike"] = strike
                    if not isinstance(euler_pole,type(None)): print("\t\t",row["comp_name"],"\n","\t\t\tEuler Pole Strike: ", strike,"\n\t\t\tPredicted Strike: ",pstrike)
                    else: print("\t\t",row["comp_name"], strike)

                #Visualize these euler poles for comparision
                if visualize:
                    pdis = np.mean(dists)
                    print("Average Distance to GC Pole: ", pdis)
                    ax = psk.plot_small_circle(plon,plat,pdis,color = "k", m=ax, geoid=Geodesic(6371.,0.), transform=ccrs.PlateCarree(), alpha=.7, linewidth=5, zorder=1)
                    if not isinstance(euler_pole,type(None)):
                        estrike = np.mean(estrikes)
                        print("Average Azimuth of Sites Relative to EP: ", estrike)
                        ep_color = plt.rcParams['axes.prop_cycle'].by_key()['color'][(ep_idx%9)+1]
                        ax = psk.plot_great_circle(euler_pole[1],euler_pole[0],estrike, m=ax, color=ep_color, geoid=Geodesic(6371.,0.), transform=ccrs.PlateCarree(), alpha=.7, linewidth=3, zorder=2)
            if visualize: #Finish of the visuals
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
            if strike < 0: strike += 360
            if strike < 180: strike += 180
            for i,row in sz_df.iterrows():
                dsk_df.at[i,"strike"] = strike
                print("\t",row["comp_name"], strike)

        else: #under determined case; just ignore
            pass

    if filter_by_quality: #If commented data was removed, add it back in
        dsk_df = dsk_df.append(bad_dsk_data)
        dsk_df.sort_values("inter_lat",inplace=True,ascending=False)

    print("--------------------------------------")
    if max_file!=None:
        #Calculate the Preturbed Poles
        (mplat,mplon,mpmag,mpa,mpb,mpphi),mpchisq,mpdof = pymax.max_likelihood_pole(pdata_splus, trial_pole=header[:3], out_path="calc_strike_splus.maxout", step=header[-1], max_steps=100, comment=comment, solve_anom_skew=False)
        (mslat,mslon,msmag,msa,msb,msphi),mschisq,msdof = pymax.max_likelihood_pole(pdata_sminus, trial_pole=header[:3], out_path="calc_strike_sminus.maxout", step=header[-1], max_steps=100, comment=comment, solve_anom_skew=False)
        print("Pole with Strikes preturbed clockwise (positive): ", (mplat,mplon,mpmag,mpa,mpb,mpphi),mpchisq,mpdof)
        print("Pole with Strikes preturbed counter-clockwise (negative): ",(mslat,mslon,msmag,msa,msb,msphi),mschisq,msdof)
        #Get their covariance matrices
        pos_geodict = geoid.Inverse(mlat,mlon,mplat,mplon)
        neg_geodict = geoid.Inverse(mlat,mlon,mslat,mslon)
        sa = max([pos_geodict["a12"],neg_geodict["a12"]])
        if pos_geodict["a12"]>neg_geodict["a12"]: max_geodict,min_geodict = pos_geodict,neg_geodict
        else: max_geodict,min_geodict = neg_geodict,pos_geodict
        sphi = max_geodict["azi1"]
        sb = min_geodict["a12"]*abs(np.sin(np.deg2rad(sphi-min_geodict["azi1"])))
        (sx,sy,sz),scov = latlon2cart(mlat,mlon,ellipse_to_cov(mlat,mlon,sa,sb,sphi))
        print("Preturbed Distances and Azimuths: ", pos_geodict["a12"],pos_geodict["azi1"],neg_geodict["a12"],neg_geodict["azi1"])
        print("Deviation of Azimuths from Anti-parallel: ", 180+(pos_geodict["azi1"]-neg_geodict["azi1"]))
        print("Strike Ellipse: ", sa, sb, sphi)
        #Add and then convert back to an ellipse
        (mlat,mlon),totcov = cart2latlon(mx,my,mz,mcov+scov)
        full_unc = cov_to_ellipse(mlat,mlon,totcov)
        print("Strike Covariance Matrix:\n",scov)
        print("Full Uncertainty: ",full_unc)

    #More Plotting stuff for comparision with Euler poles
    if not isinstance(euler_pole,type(None)):
        if visualize:
            all_strike_diffs = []
            fig_all = plt.figure(dpi=100)
            ax_all = fig_all.add_subplot(111)
            for ep_idx in range(len(strike_diffs)):
                ep_color = plt.rcParams['axes.prop_cycle'].by_key()['color'][(ep_idx%9)+1]
                #Do histogram for each individual euler pole
                print("For EP %d -> Mean, Median, Min, Max Strike Differences: "%ep_idx,sum(strike_diffs[ep_idx])/len(strike_diffs[ep_idx]),np.median(strike_diffs[ep_idx]),min(strike_diffs[ep_idx]),max(strike_diffs[ep_idx]))
                fig = plt.figure(dpi=100)
                ax = fig.add_subplot(111)
                ax.hist(strike_diffs[ep_idx],bins=np.arange(0.,4.2,0.2),color=ep_color)
                ax.axvline(sum(strike_diffs[ep_idx])/len(strike_diffs[ep_idx]),color="tab:blue",linestyle="--")
                ax.axvline(np.median(strike_diffs[ep_idx]),color="cyan")
                vis_outpath = os.path.join(os.path.dirname(dsk_path),"strike_fit_epstats_%d.png"%ep_idx)
                print("Saving: %s"%vis_outpath)
                fig.savefig(vis_outpath)
                #Do stacked histogram for all euler poles
                all_strike_diffs += list(strike_diffs[ep_idx])
                ax_all.hist(all_strike_diffs,bins=np.arange(0.,4.2,0.2),color=ep_color,zorder=len(strike_diffs)-ep_idx)
#            all_strike_diffs = reduce(lambda x,y=[]: x+y, strike_diffs)
            print("For All EP -> Mean, Median, Min, Max Strike Differences: ",sum(all_strike_diffs)/len(all_strike_diffs),np.median(all_strike_diffs),min(all_strike_diffs),max(all_strike_diffs))
            ax_all.axvline(sum(all_strike_diffs)/len(all_strike_diffs),color="tab:red",linestyle="--")
            ax_all.axvline(np.median(all_strike_diffs),color="tab:orange")
            vis_outpath = os.path.join(os.path.dirname(dsk_path),"strike_fit_all_epstats.png")
            print("Saving: %s"%vis_outpath)
            fig_all.savefig(vis_outpath)
            all_strike_diffs = reduce(lambda x,y=[]: x+y, strike_diffs)
            print("For All EP (Check) -> Mean, Median, Min, Max Strike Differences: ",sum(all_strike_diffs)/len(all_strike_diffs),np.median(all_strike_diffs),min(all_strike_diffs),max(all_strike_diffs))

    #Write out final answer in deskew file
    if isinstance(outfile,type(None)): outfile = os.path.join(os.path.dirname(dsk_path),"strike_cor_"+os.path.basename(dsk_path))
    print("Writing to %s"%str(outfile))
    utl.write_deskew_file(outfile,dsk_df)

    return full_unc

#This is to reconsile the Max file data with the Deskew file data as they may be stored different
def find_pdata_idx(comp_name, pdata, pdata_key="phs"):
    for i in range(len(pdata[pdata_key])):
        if pdata[pdata_key][i][0]==comp_name: return i

if __name__=="__main__":

    kwargs = {}
    if "-h" in sys.argv:
        help(__name__)
        sys.exit()
    if "-max" in sys.argv:
        imax = sys.argv.index("-max")
        kwargs["max_file"] = sys.argv[imax+1]
    if "-askw" in sys.argv:
        kwargs["solve_anom_skew"] = True
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


