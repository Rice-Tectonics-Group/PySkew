#!/usr/bin/env python3
"""

"""
import os,sys
import numpy as np
import pandas as pd
from scipy.stats import f
from pyskew.utilities import cov_to_ellipse
from geographiclib.geodesic import Geodesic
from functools import reduce

##############################################Main Functions################################################

def max_likelihood_pole(data, trial_pole=None, out_path=None, save_full_data_kernel=False, step=1.0, max_steps=100, convergence_level=None, comment=None, geoid=Geodesic(6371.,0.), supress_warnings=False, solve_anom_skew=False):

    #parse trial pole or generate random pole
    if isinstance(trial_pole,type(None)):
        if solve_anom_skew: trial_pole = np.array([180*np.random.rand()-90,360*np.random.rand(),1.0,0.0])
        else: trial_pole = np.array([180*np.random.rand()-90,360*np.random.rand(),1.0])
    elif ((not solve_anom_skew) and len(trial_pole)!=3): raise ValueError("Starting trial pole must be at least length 3 sequence datacontaining (lat,lon,amp). Got: %s"%str(trial_pole))
    elif (solve_anom_skew and len(trial_pole)!=4): trial_pole = [*trial_pole,0.]; print("Trial pole missing anomalous skewness adding 0 as initial guess")
    trial_pole = np.array(trial_pole) #cast as array for matrix operations

    if isinstance(convergence_level,type(None)): convergence_level=.01*step

    #decide if data is to be written to a file
    if not isinstance(out_path,type(None)):
        write_to_file = True
        fout = open(out_path,"w+")
        if comment==None: comment = "Maximum Likelyhood Pole\n"
        fout.write(comment)
        dtype_to_name = {"dec": "declinations", "inc": "inclinations",
                         "phs": "effective inclinations",
                         "ell": "poles with elliptical conf. limits",
                         "ccl": "corrected paleocolatitudes",
                         "azi": "azimuths of paleoequators",
                         "amp": "remanent amplitude factors"}
        for dtype in sorted(data.keys()):
            fout.write("the number of %s are: %d\n"%(dtype_to_name[dtype],len(data[dtype])))
        if solve_anom_skew: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment   %.2f askw\n"%tuple(trial_pole))
        else: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment\n"%tuple(trial_pole))
        fout.write("the step size is: %.2f\n"%step)
        fout.write("the convergence limit is: %.2f\n"%convergence_level)
    else: write_to_file=False

    if len(data["ell"])>0 and len(data["ell"][0][1])==5: #THIS MAY CAUSE A MUTATION PROBLEM I NEED TO LOOK INTO
        ell_1D_data = []
        for name,(plat,plon,a,b,phi) in data["ell"]:
            geodict = geoid.ArcDirect(plat,plon,phi,90)
            geodict2 = geoid.ArcDirect(plat,plon,phi+90,90)
            ell_1D_data.append([name,[[90,a/np.sqrt(2),geodict["lat2"],geodict["lon2"],geodict["azi2"]],[90,b/np.sqrt(2),geodict2["lat2"],geodict2["lon2"],geodict["azi2"]]]])
        data["ell"] = ell_1D_data

    if len(data["amp"])==0: data["amp"].append(["Fake Amplitude",[1.0,0.1,0.0,180.0,90.0]]) #if there is no amplitude data insert a fake amplitude for the inversion to avoid singular matrix

    for loop in range(max_steps):
        #empty variables for next iteration
        if solve_anom_skew: a,db = [[],[],[],[]],[]
        else: a,db = [[],[],[]],[]

        if write_to_file: fout.write("\nIteration %d\n"%loop)

        #set up the linearized eqns. for decs
        for name,(datum,sd,slat,slon) in data["dec"]:
            db.append(wrap_180(datum-linearized_declination_func(1,*trial_pole[:2],slat,slon))/sd)
            a[0].append(linearized_declination_func(2,*trial_pole[:2],slat,slon)/sd)
            a[1].append(linearized_declination_func(3,*trial_pole[:2],slat,slon)/sd)
            a[2].append(0)
            if solve_anom_skew: a[3].append(0)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for incs
        for name,(datum,sd,slat,slon) in data["inc"]:
            db.append((datum-linearized_inclination_func(1,*trial_pole[:2],slat,slon))/sd)
            a[0].append(linearized_inclination_func(2,*trial_pole[:2],slat,slon)/sd)
            a[1].append(linearized_inclination_func(3,*trial_pole[:2],slat,slon)/sd)
            a[2].append(0)
            if solve_anom_skew: a[3].append(0)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for anomaly phase shifts
        for name,(datum,sd,slat,slon,azi) in data["phs"]:
            if solve_anom_skew: db.append(wrap_180(datum-linearized_effective_inclination_func(1,*trial_pole[:2],slat,slon,azi,anom_skew=trial_pole[3]))/sd)
            else: db.append(wrap_180(datum-linearized_effective_inclination_func(1,*trial_pole[:2],slat,slon,azi))/sd)
            a[0].append(linearized_effective_inclination_func(2,*trial_pole[:2],slat,slon,azi)/sd)
            a[1].append(linearized_effective_inclination_func(3,*trial_pole[:2],slat,slon,azi)/sd)
            a[2].append(0)
            if solve_anom_skew: a[3].append(-1)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for elliptical confidence data and unoriented vertical cores
        for name,(dofs) in data["ell"]:
            for datum,sd,slat,slon,azi in dofs: #ellipse has 2 degrees of freedom iterate over both
                tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
                db.append((datum-np.rad2deg(tcol))/sd)
                a[0].append(-dot_and_deriv(2,*trial_pole[:2],slat,slon)/(np.sin(tcol)*sd))
                a[1].append(-dot_and_deriv(3,*trial_pole[:2],slat,slon)/(np.sin(tcol)*sd))
                a[2].append(0)
                if solve_anom_skew: a[3].append(0)
                if write_to_file and save_full_data_kernel:
                    if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                    a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                    fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for corrected colatitudes from unoriented vertical cores
        for name,(datum,sd,slat,slon) in data["ccl"]:
            tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
            db.append((datum-np.rad2deg(tcol))/sd)
            a[0].append(-dot_and_deriv(2,*trial_pole[:2],slat,slon)/(np.sin(tcol)*sd))
            a[1].append(-dot_and_deriv(3,*trial_pole[:2],slat,slon)/(np.sin(tcol)*sd))
            a[2].append(0)
            if solve_anom_skew: a[3].append(0)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for azimuths of paleoequators
        for name,(datum,sd,slat,slon) in data["azi"]:
            db.append((datum-linearized_declination_func(1,*trial_pole[:2],slat,slon))/sd)
            a[0].append(linearized_declination_func(2,*trial_pole[:2],slat,slon)/sd)
            a[1].append(linearized_declination_func(3,*trial_pole[:2],slat,slon)/sd)
            a[2].append(0)
            if solve_anom_skew: a[3].append(0)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #set up the linearized eqns. for remanent amplitude factors
        for name,(datum,sd,slat,slon,azi) in data["amp"]:
            db.append((datum-linearized_amplitude_func(1,*trial_pole[:3],slat,slon,azi))/sd)
            a[0].append(linearized_amplitude_func(2,*trial_pole[:3],slat,slon,azi)/sd)
            a[1].append(linearized_amplitude_func(3,*trial_pole[:3],slat,slon,azi)/sd)
            a[2].append(linearized_amplitude_func(4,*trial_pole[:3],slat,slon,azi)/sd)
            if solve_anom_skew: a[3].append(0)
            if write_to_file and save_full_data_kernel:
                if len(a[0])-1==0: fout.write("Data Kernal and normalized residual:\td/dplat\td/dplon\td/dm\n")
                a_str = reduce(lambda x,y: x+"\t"+y, ["%.2f"%a[i][-1] for i in range(len(a))])
                fout.write("                             row %d:	%s\t%.4f\n"%(len(a[0])-1,a_str,db[-1]))

        #Solve the linearized equations (a.T @ a) @ dm = a.T @ db
        a,db = np.array(a).T,np.array(db)

        #calculate chi-square for this trial
        chisq = db.T @ db

        #calculate model update
        try: cov = np.linalg.inv(a.T @ a)
        except np.linalg.LinAlgError:
            if not supress_warnings:
                print("******************************************")
                print("Singular Matrix")
                print(data, trial_pole)
                print(a)
                print("******************************************")
            cov = np.zeros([3,3])
            break
        dm = ((cov @ a.T) @ db)*step

        if write_to_file: #Write out this iteration
            if solve_anom_skew: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment   %.2f askw\n"%tuple(trial_pole))
            else: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment\n"%tuple(trial_pole))
            fout.write("Covariance Matrix:\n")
            if solve_anom_skew:
                fout.write("\tlat\tlon\tm\taskw\nlat\t%.4f\t%.4f\t%.4f\t%.4f\nlon\t%.4f\t%.4f\t%.4f\t%.4f\nm\t%.4f\t%.4f\t%.4f\t%.4f\naskw\t%.4f\t%.4f\t%.4f\t%.4f\n"%tuple(cov.flatten().astype(float)))
                fout.write("Model Step: %.2f N   %.2f E  %.2f mag. moment   %.2f askw\n"%tuple(dm))
            else:
                fout.write("\tlat\tlon\tm\nlat\t%.4f\t%.4f\t%.4f\nlon\t%.4f\t%.4f\t%.4f\nm\t%.4f\t%.4f\t%.4f\n"%tuple(cov.flatten().astype(float)))
                fout.write("Model Step: %.2f N   %.2f E  %.2f mag. moment\n"%tuple(dm))
            fout.write("chi squared = %.2f\n"%chisq)

        #update the model
        trial_pole += dm
        #handle roll over error in lat and lon
        if abs(trial_pole[0])>90:
            if write_to_file:
                fout.write("******************************************\n")
                fout.write("the model has crossed over a pole\n")
                fout.write("******************************************\n")
            while abs(trial_pole[0])>90:
                trial_pole[0] = 180*np.sign(trial_pole[0])-trial_pole[0]
                trial_pole[1] = -180*np.sign(trial_pole[1])+trial_pole[1]
        if trial_pole[1]>360: trial_pole[1] = trial_pole[1]%360
        elif trial_pole[1]<-360: trial_pole[1] = trial_pole[1]%-360

        #handle fixed anomalous skewness case
        if ((isinstance(solve_anom_skew,int) and not isinstance(solve_anom_skew,bool)) or isinstance(solve_anom_skew,float)) and solve_anom_skew>1:
            trial_pole[-1] = solve_anom_skew
            dm[-1] = 0.

        #test convergence criteria
        if (sum(abs(dm[:-1]))+abs(.1*dm[-1]))<convergence_level:
            break

    if not supress_warnings and loop==max_steps-1:
        print("******************************************")
        print("Max did not converge after %d iterations"%max_steps)
#        print(data,trial_pole)
#        print(dm)
#        print(a)
#        print(db)
#        print((sum(abs(dm[:-1]))+abs(.1*dm[-1])),convergence_level)
        print("******************************************")

    #figure out uncertainties
    maj_se,min_se,phi = cov_to_ellipse(*trial_pole[:2],cov)
    if solve_anom_skew: askw_se = np.sqrt(cov[3,3])
    dof = len(data["dec"])+len(data["inc"])+len(data["phs"])+2*len(data["ell"])+len(data["ccl"])+len(data["azi"])+len(data["amp"])-len(trial_pole)

    if write_to_file:
        fout.write("\n\n\nTotal Number of Iterations: %d\n\n"%(loop+1))
        if solve_anom_skew: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment   %.2f askw\n"%tuple(trial_pole))
        else: fout.write("the trial pole is: %.2f N   %.2f E  %.2f mag. moment\n"%tuple(trial_pole))
        f_factor = f.ppf(.95,2,dof)
        if solve_anom_skew: fout.write("The 95%" + " confidence ellipse is:\n\tMajor Semi-Axis = %.2fº\n\tMinor Semi-Axis = %.2fº\n\tAzimuth of Major Axis = %.2f\n\tAnomalous Skewness Unc. = %.3f\n"%(np.sqrt(f_factor)*maj_se,np.sqrt(f_factor)*min_se,phi,np.sqrt(f_factor)*askw_se))
        else: fout.write("The 95%" + " confidence ellipse is:\n\tMajor Semi-Axis = %.2fº\n\tMinor Semi-Axis = %.2fº\n\tAzimuth of Major Axis = %.2f\n"%(np.sqrt(f_factor)*maj_se,np.sqrt(f_factor)*min_se,phi))
        fout.write("The F-ratio factor = %.4f\n\n"%f_factor)
        fout.write("chi squared = %.2f with %d degrees of freedom\n"%(chisq,dof))
        fout.write("reduced chi squared = %.2f\n\n"%(chisq/dof))
        importances = sum(cov[i,i]*a[:,i]**2 if i==j else cov[i,j]*a[:,i]*a[:,j] for i in range(3) for j in range(3))
        write_out_data_table(fout,data,trial_pole,cov,importances,geoid=geoid)

    return (*trial_pole,maj_se,min_se,phi),chisq,dof

def get_chisq_no_magnitude(data,trial_pole=[0.0,0.0,1.0],geoid=Geodesic(6371000.,0.)):
    trial_pole = np.array(trial_pole) #cast as array for matrix operations

    if len(data["ell"])>0 and len(data["ell"][0][1])==5: #THIS MAY CAUSE A MUTATION PROBLEM I NEED TO LOOK INTO
        ell_1D_data = []
        for name,(plat,plon,a,b,phi) in data["ell"]:
            geodict = geoid.ArcDirect(plat,plon,phi,90)
            geodict2 = geoid.ArcDirect(plat,plon,phi+90,90)
            ell_1D_data.append([name,[[90,a/np.sqrt(2),geodict["lat2"],geodict["lon2"],geodict["azi2"]],[90,b/np.sqrt(2),geodict2["lat2"],geodict2["lon2"],geodict["azi2"]]]])
        data["ell"] = ell_1D_data

    #empty variables for next iteration
    db = []

    #set up the linearized eqns. for decs
    for name,(datum,sd,slat,slon) in data["dec"]:
        db.append(wrap_180(datum-linearized_declination_func(1,*trial_pole[:2],slat,slon))/sd)

    #set up the linearized eqns. for incs
    for name,(datum,sd,slat,slon) in data["inc"]:
        db.append((datum-linearized_inclination_func(1,*trial_pole[:2],slat,slon))/sd)

    #set up the linearized eqns. for anomaly phase shifts
    for name,(datum,sd,slat,slon,azi) in data["phs"]:
        db.append(wrap_180(datum-linearized_effective_inclination_func(1,*trial_pole[:2],slat,slon,azi))/sd)

    #set up the linearized eqns. for elliptical confidence data and unoriented vertical cores
    for name,(dofs) in data["ell"]:
        for datum,sd,slat,slon,azi in dofs: #ellipse has 2 degrees of freedom iterate over both
            tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
            db.append((datum-np.rad2deg(tcol))/sd)

    #set up the linearized eqns. for corrected colatitudes from unoriented vertical cores
    for name,(datum,sd,slat,slon) in data["ccl"]:
        tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
        db.append((datum-np.rad2deg(tcol))/sd)

    #set up the linearized eqns. for azimuths of paleoequators
    for name,(datum,sd,slat,slon) in data["azi"]:
        db.append((datum-linearized_declination_func(1,*trial_pole[:2],slat,slon))/sd)

    #calculate chi-square for this trial
    db = np.array(db)
    return db.T @ db

def chisq_surface(data,step=1.0,lat_range=(-90,90),lon_range=(0,360),geoid=Geodesic(6371000.,0.)):
    lats = np.arange(lat_range[0],lat_range[1]+step,step)
    lons = np.arange(lon_range[0],lon_range[1]+step,step)

    chi2_grid = np.zeros([len(lats),len(lons)])
    for i,lat in enumerate(lats):
        for j,lon in enumerate(lons):
            chi2_grid[i,j] = get_chisq_no_magnitude(data, trial_pole=[lat,lon,1.0])

    return lons,lats,chi2_grid

##############################################I/O Functions################################################

def read_max_file(max_path,sep=","):

    with open(max_path,'r') as fin:

        #read header
        comment = fin.readline()
        header = list(map(float,fin.readline().split(sep)))
        tplat,tplon,tpm,ndir,ninc,naps,nell,ncorcl,naz,nraf,step = header

        #build data
        data = {"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]}

        #read dec/inc pair
        for i in range(int(ndir)): #expects two lines with (datum, standard err., site lat, site lon) dec first
            name_datum = fin.readline().rstrip("\n")
            data["dec"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])
            name_datum = fin.readline().rstrip("\n")
            data["inc"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

        #read inc only
        for i in range(int(ninc)): #expects one line with (datum, standard err., site lat, site lon)
            name_datum = fin.readline().rstrip("\n")
            data["inc"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

        #read anomaly phase shift data
        for i in range(int(naps)): #expects one line with (datum, standard err., site lat, site lon, azimuth of isochron)
            name_datum = fin.readline().rstrip("\n")
            data["phs"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

        #read elliptical confidence data
        for i in range(int(nell)): #expects one line with (plat, plon, major semi-axis, minor semi-axis, azimuth of major axis)
            name_datum = fin.readline().rstrip("\n")
            data["ell"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

        #read corrected colatitudes from unoriented vertical cores
        #(assumes they've been corrected already according to Cox & Gordon 1984)
        for i in range(int(ncorcl)): #expects one line with (datum, standard err., site lat, site lon)
            name_datum = fin.readline().rstrip("\n")
            data["ccl"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

        #read in the azimuths of paleoequators
        for i in range(int(naz)): #expects one line with (datum, standard err., site lat, site lon)
            name_datum = fin.readline().rstrip("\n")
            data["azi"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])
            data["azi"][1][0] -= 90

        #read in remanent amplitude factors
        for i in range(int(nraf)): #expects one line with (datum, standard err., site lat, site lon, azi)
            name_datum = fin.readline().rstrip("\n")
            data["amp"].append([name_datum,[float(item) for item in fin.readline().split(sep)]])

    return comment,header,data

def read_csv(csv_path,sep="\t",geoid=Geodesic(6371000.,0.)):
    df = pd.read_csv(csv_path,sep=sep)
    return parse_apwp_df_or_list_to_max(df,geoid=geoid)

def to_csv(out_path,sep="\t"):
    pass

def parse_apwp_df_or_list_to_max(df_or_lst,geoid=Geodesic(6371000.,0.)):
    if isinstance(df_or_lst,list): iterable = df_or_lst
    else: iterable = df_or_lst.to_dict("record")

    data = {"dec":[],"inc":[],"phs":[],"ell":[],"ccl":[],"azi":[],"amp":[]}

    for row in iterable:
        if "DecOnly" in row.keys() and row["DecOnly"]:
            plon,plat,slon,slat,error = row["Longitude"],row["Latitude"],row["Site Longitude"],row["Site Latitude"],row["Major (SE)"]
            geodict = geoid.Inverse(slat,slon,plat,plon)
            azi_mean = geodict["azi1"]
            geodict = geoid.ArcDirect(plat,plon,geodict['azi2']-90,error)
            geodict = geoid.Inverse(slat,slon,geodict['lat2'],geodict['lon2'])
            azi_se = (azi_mean-geodict['azi1'])/np.sqrt(2)
            data["dec"].append([str(row["Name"])+" "+str(row["Age"]),[azi_mean,azi_se,slat,slon]])
        elif "IncOnly" in row.keys() and row["IncOnly"]:
            if "PaleoCoLatitude" in row.keys() and not np.isnan(row["PaleoCoLatitude"]): datum,se = row["PaleoCoLatitude"],row["Major (SE)"]
            else:
                geodict = geoid.Inverse(row["Latitude"],row["Longitude"],row["Site Latitude"],row["Site Longitude"])
                datum = geodict["a12"]
                a,b,phi,azi = row["Major (SE)"],row["Minor (SE)"],row["Azimuth"],geodict["azi1"]
                se = ((a*b)/np.sqrt((np.cos(np.deg2rad(phi-azi))*b)**2 + (np.sin(np.deg2rad(phi-azi))*a)**2))/np.sqrt(2)
            data["ccl"].append([str(row["Name"])+" "+str(row["Age"]),[datum,se,row["Site Latitude"],row["Site Longitude"]]])
        else:
            data["ell"].append([str(row["Name"])+" "+str(row["Age"]),[row["Latitude"],row["Longitude"],row["Major (SE)"],row["Minor (SE)"],row["Azimuth"]]])

    return data

def write_out_data_table(fout,data,trial_pole,cov,importances,geoid=Geodesic(6371000.,0.)):

    it = 0

    fout.write("Declinations:\n")
    fout.write("name\tlat\tlon\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon) in data["dec"]:
        model_datum = linearized_declination_func(1,*trial_pole[:2],slat,slon)
        residual = wrap_180(datum-model_datum)
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,datum,sd,model_datum,residual,importances[it]))
        it+=1

    fout.write("Inclinations:\n")
    fout.write("name\tlat\tlon\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon) in data["inc"]:
        model_datum = linearized_inclination_func(1,*trial_pole[:2],slat,slon)
        residual = datum-model_datum
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,datum,sd,model_datum,residual,importances[it]))
        it+=1

    fout.write("Effective Inclinations:\n")
    fout.write("name\tlat\tlon\tazi\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon,azi) in data["phs"]:
        model_datum = linearized_effective_inclination_func(1,*trial_pole[:2],slat,slon,azi)
        residual = datum-model_datum
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,azi,datum,sd,model_datum,residual,importances[it]))
        it+=1

    fout.write("Poles:\n")
    fout.write("name\taxis\tlat\tlon\tstd err\tazi\tmodel lat\tmodel lon\tresidual\timportance\n")
    for name,(dofs) in data["ell"]:
        #Major Axis
        datum,sd,slat,slon,azi = dofs[0]
        tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
        maj_resid = (datum-np.rad2deg(tcol))
        major = np.sqrt(2)*sd
        geodict = geoid.ArcDirect(slat,slon,azi,datum)
        fout.write("%s\tmajor\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,geodict["lat2"],geodict["lon2"],major,geodict["azi2"],trial_pole[0],trial_pole[1],maj_resid,importances[it]))
        it+=1
        #Minor Axis
        datum,sd,slat,slon,azi = dofs[0]
        tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
        min_resid = (datum-np.rad2deg(tcol))
        minor = np.sqrt(2)*sd
        geodict = geoid.ArcDirect(slat,slon,azi,datum)
        fout.write("\tminor\t\t\t%.2f\t\t\t\t%.2f\t%.3f\n"%(minor,min_resid,importances[it]))
        it+=1

    fout.write("Paleo Colatitudes:\n")
    fout.write("name\tlat\tlon\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon) in data["ccl"]:
        tcol = np.arccos(dot_and_deriv(1,*trial_pole[:2],slat,slon))
        model_datum = np.rad2deg(tcol)
        residual = datum-model_datum
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,datum,sd,model_datum,residual,importances[it]))
        it+=1

    fout.write("Azimuths of Paleoequators:\n")
    fout.write("name\tlat\tlon\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon) in data["azi"]:
        model_datum = linearized_declination_func(1,*trial_pole[:2],slat,slon)
        residual = datum-model_datum
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,datum,sd,model_datum,residual,importances[it]))
        it+=1

    fout.write("Remanent Amplitude Factors:\n")
    fout.write("name\tlat\tlon\tazi\tdatum\tstd err\tmodel\tresidual\timportance\n")
    for name,(datum,sd,slat,slon,azi) in data["amp"]:
        model_datum = linearized_amplitude_func(1,*trial_pole[:3],slat,slon,azi)
        residual = datum-model_datum
        fout.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n"%(name,slat,slon,azi,datum,sd,model_datum,residual,importances[it]))
        it+=1

##############################################Utility Functions################################################

def dot_and_deriv(idx,plat,plon,slat,slon):
    rplat,rplon,rslat,rslon = np.deg2rad([plat,plon,slat,slon])
    if idx==1: #dot product of pole and site
        return np.cos(rplat)*np.cos(rplon)*np.cos(rslat)*np.cos(rslon) + np.cos(rplat)*np.sin(rplon)*np.cos(rslat)*np.sin(rslon) + np.sin(rplat)*np.sin(rslat)
    elif idx==2: #d(p @ s)/d(plat)
        return -np.sin(rplat)*np.cos(rplon)*np.cos(rslat)*np.cos(rslon) - np.sin(rplat)*np.sin(rplon)*np.cos(rslat)*np.sin(rslon)+np.cos(rplat)*np.sin(rslat)
    elif idx==3: #d(p @ s)/d(plon)
        return -np.cos(rplat)*np.sin(rplon)*np.cos(rslat)*np.cos(rslon) + np.cos(rplat)*np.cos(rplon)*np.cos(rslat)*np.sin(rslon)
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def linearized_declination_func(idx,plat,plon,slat,slon):
    rplat,rplon,rslat,rslon = np.deg2rad([plat,plon,slat,slon])
    if idx==1: #declination given spherical earth approx. in degrees
        return (360-np.rad2deg(np.arctan2(np.cos(rplat)*np.sin(rslon-rplon),-np.sin(rslat)*np.cos(rplat)*np.cos(rslon-rplon)+np.cos(rslat)*np.sin(rplat))))%360
    elif idx==2: #d(dec)/d(plat)
        tmp = 1-dot_and_deriv(1,plat,plon,slat,slon)**2
        return np.cos(rslat)*np.sin(rslon-rplon)/tmp #order of opperations may be wrong idk if it's different between python and fortran
    elif idx==3: #d(dec)/d(plon)
        tmp = 1-dot_and_deriv(1,plat,plon,slat,slon)**2
        return -(np.cos(rplat)*np.sin(rslat)-np.sin(rplat)*np.cos(rslat)*np.cos(rslon-rplon))*np.cos(rplat)/tmp
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def linearized_inclination_func(idx,plat,plon,slat,slon):
    gg = dot_and_deriv(1,plat,plon,slat,slon)
    if idx==1: #inclination given spherical earth approx. in degrees
        return np.rad2deg(np.arctan2(2*gg,np.sqrt(1-gg**2)))
    elif idx==2: #d(inc)/d(plat)
        return (2*dot_and_deriv(2,plat,plon,slat,slon))/((1+3*(gg**2))*np.sqrt(1-gg**2))
    elif idx==3: #d(inc)/d(plon)
        return (2*dot_and_deriv(3,plat,plon,slat,slon))/((1+3*(gg**2))*np.sqrt(1-gg**2))
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def linearized_effective_inclination_func(idx,plat,plon,slat,slon,azi,anom_skew=0):
    gg = dot_and_deriv(1,plat,plon,slat,slon)
    d_azi = np.deg2rad(azi-linearized_declination_func(1,plat,plon,slat,slon))
    if idx==1: #effective inclination given spherical earth approx. in degrees
        return np.rad2deg(np.arctan2(2*gg,np.sin(d_azi)*np.sqrt(1-gg**2))) - anom_skew
    elif idx==2: #d(eff_inc)/d(plat)
        dd = linearized_declination_func(2,plat,plon,slat,slon)
        dg = dot_and_deriv(2,plat,plon,slat,slon)
        return (2/(np.sin(d_azi)**2+(4-np.sin(d_azi)**2)*(gg**2)))*((np.sin(d_azi)*dg)/np.sqrt(1-gg**2) + gg*np.sqrt(1-gg**2)*np.cos(d_azi)*dd)
    elif idx==3: #d(eff_inc)/d(plon)
        dd = linearized_declination_func(3,plat,plon,slat,slon)
        dg = dot_and_deriv(3,plat,plon,slat,slon)
        return (2/(np.sin(d_azi)**2+(4-np.sin(d_azi)**2)*(gg**2)))*((np.sin(d_azi)*dg)/np.sqrt(1-gg**2) + gg*np.sqrt(1-gg**2)*np.cos(d_azi)*dd)
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def linearized_amplitude_func(idx,plat,plon,pm,slat,slon,azi):
    daf = dipole_amplitude_factor(1,plat,plon,slat,slon)
    raf = remanent_amplitude_factor(1,plat,plon,slat,slon,azi)
    if idx==1: #model total amplitude factor
        return pm*daf*raf
    elif idx==2 or idx==3: #d(af)/d(plat) or d(af)/d(plon)
        return pm*(daf*remanent_amplitude_factor(idx,plat,plon,slat,slon,azi)+raf*dipole_amplitude_factor(idx,plat,plon,slat,slon))
    elif idx==4: # d(af)/d(pm)
        return daf*raf
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def dipole_amplitude_factor(idx,plat,plon,slat,slon):
    gg = dot_and_deriv(1,plat,plon,slat,slon)
    if idx==1: #dipole amplitude factor
        return np.sqrt(1+3*gg**2)
    elif idx==2 or idx==3: #d(daf)/d(plat) or d(daf)/d(plon)
        return np.deg2rad(((3*gg)/np.sqrt(1+3*gg**2))*dot_and_deriv(idx,plat,plon,slat,slon))
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def remanent_amplitude_factor(idx,plat,plon,slat,slon,azi):
    inc = np.deg2rad(linearized_inclination_func(1,plat,plon,slat,slon))
    rei = np.deg2rad(linearized_effective_inclination_func(1,plat,plon,slat,slon,azi)) #remanent effective inclination
    if idx==1: #remanent amplitude factor
        return np.sin(inc)/np.sin(rei)
    elif idx==2 or idx==3: #d(raf)/d(plat) or d(raf)/d(plon)
        tmp1 = np.cos(inc)*linearized_inclination_func(idx,plat,plon,slat,slon)/np.sin(rei)
        tmp2 = np.sin(inc)*np.cos(rei)*linearized_effective_inclination_func(idx,plat,plon,slat,slon,azi)/(np.sin(rei)**2)
        return np.deg2rad(tmp1-tmp2)
    else: raise ValueError("No case %d in function first value must be on interval [1:3]"%idx)

def wrap_180(val):
    if abs(val)>=180: return val-360*np.sign(val)
    else: return val

def wrap_360(val):
    return (360+val)%360


#Scripting portion of the code
if __name__=="__main__":
    if "-h" in sys.argv or "--help" in sys.argv: help(__name__)

    if len(sys.argv)<3: raise ValueError("Requires a max input file and output file. For formating information run with -h")

    if "-max_iter" in sys.argv:
        idx = sys.argv.index("-max_iter")
        try: max_iter = int(sys.argv[idx+1])
        except ValueError: raise ValueError("Maximum iterations must be a integer. Recieved: %s"%sys.argv[idx+1])
    else: max_iter = 100

    if "-v" in sys.argv: verbose = True
    else: verbose = False

    if "-askw" in sys.argv: solve_anom_skew = True
    else: solve_anom_skew = False

    comment,header,data = read_max_file(sys.argv[1])
    max_likelihood_pole(data, trial_pole=header[:3],out_path=sys.argv[2], save_full_data_kernel=verbose,step=header[-1],max_steps=max_iter,comment=comment, solve_anom_skew=solve_anom_skew)




