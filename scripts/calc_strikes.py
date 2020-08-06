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

def calc_strikes_and_add_err(dsk_path,mlat=90,mlon=0,ma=1,mb=1,mphi=0,geoid=Geodesic(6371,0.0),outfile=None,filter_by_quality=False):
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
                data["ccl"].append([row["comp_name"],[90.0,0.10,row["inter_lat"],row["inter_lon"]]])
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data)
            for i in range(len(data["ccl"])):
                data["ccl"][i][1][1] *= np.sqrt(chisq)
            (plat,plon,_,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data)
            print("\t",(plat,plon,maj_se,min_se,phi),chisq,dof)
            scov = ellipse_to_cov(plat,plon,maj_se,min_se,phi)
            tcov += scov
            for i,row in sz_df.iterrows():
                strike = geoid.Inverse(plat,plon,row["inter_lat"],row["inter_lon"])["azi2"]+90
                if strike < 0: strike += 180
                dsk_df.at[i,"strike"] = strike
                print("\t\t",row["comp_name"], strike)
        else: #under or equal determined case
            strike = geoid.Inverse(sz_df.iloc[0]["inter_lat"],sz_df.iloc[0]["inter_lon"],sz_df.iloc[1]["inter_lat"],sz_df.iloc[1]["inter_lon"])["azi1"]
            if strike < 0: strike += 180
            for i,row in sz_df.iterrows():
                dsk_df.at[i,"strike"] = strike
                print("\t",row["comp_name"], strike)
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

    calc_strikes_and_add_err(sys.argv[1],**kwargs)
