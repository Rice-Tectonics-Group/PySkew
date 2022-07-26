import os, glob
import shutil
import subprocess
import pandas as pd
import numpy as np
import statistics as stat
import pmagpy.ipmag as ipmag
from geographiclib.geodesic import Geodesic
import pyskew.geographic_preprocessing as pg
import pyskew.utilities as utl
from scipy.signal import hilbert, butter, filtfilt

def filter_deskew_and_calc_aei(deskew_path, spreading_rate_path=None, anomalous_skewness_model_path=None):
    """Creates Datatable"""

    deskew_df = utl.open_deskew_file(deskew_path)
    asf,srf,sz_list = get_asf_srf(spreading_rate_path,anomalous_skewness_model_path)

    return calc_aei(deskew_df,srf,asf)

def calc_aei(deskew_df,srf,asf):

    deskew_df["ei"] = [(90. if ".Vd." in comp else 0.) for ps,comp in zip(deskew_df["phase_shift"],deskew_df["comp_name"])]

    deskew_df["aei"] = [180. - row['phase_shift']- row["ei"] + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2)) for i,row in deskew_df.iterrows()]

    for i,row in deskew_df[deskew_df['track_type']=='ship'].iterrows():
        decimal_year = get_shipmag_decimal_year(row)
        if decimal_year==None: print("Intersection point could not be found in data file so IGRF could not be calculated and aei could not be found please check your data, skipping %s"%row['comp_name']); continue
        igrf = ipmag.igrf([decimal_year,0,float(row['inter_lat']),float(row['inter_lon'])])
        alpha = float(row['strike']) - igrf[0]
        e = np.rad2deg(np.arctan2(np.tan(np.deg2rad(igrf[1])),np.sin(np.deg2rad(alpha))))
        aei = 180 - e - float(row['phase_shift']) + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2))
        deskew_df.at[i,'ei'] = e
        deskew_df.at[i,'aei'] = aei
    return deskew_df

def row_calc_aei(row,srf,asf):
    if row["track_type"] == "aero":
        row["ei"] = 90. if ".Vd." in row["comp_name"] else 0.
        row["aei"] = 180. - float(row['phase_shift']) - row["ei"] + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2))
    else:
        decimal_year = get_shipmag_decimal_year(row)
        if decimal_year==None: raise ValueError("Intersection point could not be found in data file so IGRF could not be calculated and aei could not be found please check your data, skipping %s"%row['comp_name'])
        igrf = ipmag.igrf([decimal_year,0,float(row['inter_lat']),float(row['inter_lon'])])
        alpha = float(row['strike']) - igrf[0]
        e = np.rad2deg(np.arctan2(np.tan(np.deg2rad(igrf[1])),np.sin(np.deg2rad(alpha))))
        aei = 180. - e - float(row['phase_shift']) + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2))
        row['ei'] = e
        row['aei'] = aei
    return row

def get_asf_srf(spreading_rate_path=None,anomalous_skewness_model_path=None):
    if spreading_rate_path==None:
        if os.path.isfile(os.path.join('..','raw_data','spreading_rate_model.txt')):
            spreading_rate_path = os.path.join('..','raw_data','spreading_rate_model.txt')
        else:
            spreading_rate_path = input("spreading rate model file could not be found in raw_data/spreading_rate_model.txt please provide a path to a spreading rate model file here: ")
    if anomalous_skewness_model_path==None:
        if os.path.isfile(os.path.join('..','raw_data','anomalous_skewness_model.txt')):
            anomalous_skewness_model_path = os.path.join('..','raw_data','anomalous_skewness_model.txt')
        else:
            anomalous_skewness_model_path = input("anomalous skewness model file could not be found in raw_data/anomalous_skewness_model.txt please provide a path to a anomalous skewness model file here: ")

    srf,sz_list = generate_spreading_rate_model(spreading_rate_path)
    asf = generate_anomalous_skewness_model(anomalous_skewness_model_path)

    return asf,srf,sz_list

def phase_shift_data(mag_data,phase_shift):
    mag = np.array(list(map(float,mag_data)))

    N1 = len(mag)

    #   Add a head and a tail for the data
    #   so the magnetic anamaly profile will have two smooth ends

    # Linear Tappers
#    magpre=np.arange(0,1,0.001)*mag[0]
#    magpost=np.arange(.999,-.001,-.001)*mag[-1]

    #Cosine Tappers
    magpre=((np.cos(np.linspace(-np.pi,0,1000))+1)/2)*mag[0]
    magpost=((np.cos(np.linspace(0,np.pi,1000))+1)/2)*mag[-1]

    magnew=np.concatenate([magpre,mag,magpost])

    #   Fourier transform
    N2=len(magnew)
    MAG=np.fft.fft(magnew)

    #   Apply the filter
    phi=np.deg2rad(float(phase_shift))
    if N2%2==0:
        MAG2_head=MAG[0:int(N2/2)]*np.exp(-1j*phi)
        MAG2_tail=MAG[int(N2/2):N2]*np.exp(1j*phi)
        MAG2 = np.concatenate([MAG2_head,MAG2_tail])
    else:
        MAG2_head=MAG[0:int((N2-1)/2)]*np.exp(-1j*phi)
        MAG2_tail=MAG[int((N2)/2):N2]*np.exp(1j*phi)
        MAG2 = np.concatenate([MAG2_head,MAG2_tail])

    #   Inverse Fourier transform
    NewMag=np.real(np.fft.ifft(MAG2))

    #   Truncate the head and the tail and return the data
    return NewMag[1000:N1+1000]

def correct_site(site_cor_path,deskew_path,dist_e=.5):
    #backup the .deskew file
    if not os.path.isfile(deskew_path+'.ccbak'):
        print("backing up %s to %s"%(deskew_path,deskew_path+'.ccbak'))
        shutil.copyfile(deskew_path,deskew_path+'.ccbak')

    #read in the deskew and site_cor file
    deskew_df = utl.open_deskew_file(deskew_path)
    site_cor_df = pd.read_csv(site_cor_path,sep="\t")

    #copy the deskew df so I can change it and save the corrected latitudes and longitudes
    new_deskew_df = deskew_df.copy()
    for i,drow in deskew_df.iterrows():
        crow = site_cor_df[site_cor_df["comp_name"]==drow["comp_name"]]

        if drow['comp_name'].startswith('#'): continue #commented lines check
        if crow.empty: print("no correction found for component %s"%drow["comp_name"]); continue #no correction for this component check

        if drow['track_type'] == 'aero':
            #Find other component direction so we can average the shift between components
            if 'E' in drow["comp_name"]:
                other_crow = site_cor_df[site_cor_df["comp_name"]==drow["comp_name"].replace('E','V')]
            elif 'V' in drow["comp_name"]:
                other_crow = site_cor_df[site_cor_df["comp_name"]==drow["comp_name"].replace('V','E')]
            else: print("Problem determining component direction for %s"%drow["comp_name"]); continue
            #check that the intersept distance correction between E and V are not more than 3 deg different
            if abs(float(crow['correction']) - float(other_crow['correction']))>3:
                print("correction for %s is >3 km different from the other componenet's correction, and the average may be off"%(drow['comp_name']))
            correction = (float(crow['correction'])+float(other_crow['correction']))/2
        elif drow['track_type'] == 'ship':
            correction = float(crow['correction'])
        else:
            print("could not determine the track type for %s please check your deskew file, skipping"%drow["comp_name"]); continue

        corrected_lon,corrected_lat,corrected_dist = get_lon_lat_from_plot_pick(drow,correction,dist_e=dist_e)

        new_deskew_df.at[i, 'inter_lat'] =  corrected_lat
        new_deskew_df.at[i, 'inter_lon'] =  corrected_lon

    new_deskew_df.to_csv(deskew_path,sep="\t",index=False)

def old_get_lon_lat_from_plot_pick(deskew_row,plot_pick,dist_e=.01,verbose=False):
    drow,correction=deskew_row,plot_pick

    data_file_path = os.path.join(drow["data_dir"],drow["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","idk","mag","lat","lon"],delim_whitespace=True)

    projected_distances = utl.calc_projected_distance(drow['inter_lon'],drow['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),drow['strike'])

    found_dist=False
    for j,row in projected_distances.iterrows():
        if row['dist']>=correction-dist_e and row['dist']<=correction+dist_e:
            picked_lat = round(row['lat'],3) #change lat for the new deskew file
            picked_lon = round(utl.convert_to_0_360(row['lon']),3) #change lon for the new deskew file
            picked_distance = row['dist']
            found_dist=True
            break

    if found_dist:
        if verbose: print("found lat lon of %s at a distance %.3f"%(drow["comp_name"],picked_distance))
    else:
        if verbose: print("couldn't find picked distance in datafile to calculate lat and lon for %s"%drow["comp_name"])
        return (drow['inter_lon'],drow['inter_lat'],0)

    return picked_lon,picked_lat,picked_distance

def old_get_idx_from_plot_pick(deskew_row,plot_pick,dist_e=.01):
    drow,correction=deskew_row,plot_pick

    data_file_path = os.path.join(drow["data_dir"],drow["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","idk","mag","lat","lon"],delim_whitespace=True)

    projected_distances = utl.calc_projected_distance(drow['inter_lon'],drow['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),drow['strike'])

    found_dist=False
    for j,row in projected_distances.iterrows():
        if row['dist']>=correction-dist_e and row['dist']<=correction+dist_e:
            picked_idx = j
            picked_distance = row['dist']
            found_dist=True
            break

    if found_dist:
        print("found lat lon of %s at a distance %.3f"%(drow["comp_name"],picked_distance))
    else:
        print("couldn't find picked distance in datafile to calculate lat and lon for %s"%drow["comp_name"]); return (drow['inter_lon'],drow['inter_lat'],0)

    return picked_idx,picked_distance

def get_lon_lat_from_plot_pick(deskew_row,plot_pick,flip=False,dist_e=None):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = utl.open_mag_file(data_file_path)

    if flip: projected_distances = utl.calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),(180+deskew_row['strike'])%360)
    else: projected_distances = utl.calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),deskew_row['strike'])

    min_idx = (projected_distances["dist"]-plot_pick).abs().idxmin()

    return projected_distances["lon"][min_idx],projected_distances["lat"][min_idx],projected_distances["dist"][min_idx]

def get_idx_from_plot_pick(deskew_row,plot_pick,flip=False,dist_e=None):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = utl.open_mag_file(data_file_path)

    if flip: projected_distances = utl.calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),(180+deskew_row['strike'])%360)
    else: projected_distances = utl.calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),deskew_row['strike'])

    min_idx = (projected_distances["dist"]-plot_pick).abs().idxmin()

    return min_idx,projected_distances["dist"][min_idx]

def create_maxtab_file(deskew_path,anomoly_name,outfile=None):
    deskew_df = utl.open_deskew_file(deskew_path)
    dates_data = pd.read_csv("../raw_data/dates.aeromag",sep='\t',header=0,index_col=0)
    out_str = ""
    for i,row in deskew_df.iterrows():
        if '#' in row['comp_name']: continue
        track_name = row['comp_name'].split('.')[0]
        if row['track_type']=='ship': date = "%.2f"%get_shipmag_decimal_year(row)
        else: date = "%.2f"%float(dates_data.loc[track_name]["decimal_year"])
        phase_shift = "%.2f"%utl.convert_to_0_360(float(row['phase_shift']))
        out_str += ' '*(17-len(row['comp_name']))+ row['comp_name'] + ' '
        out_str += 'V' if 'Vd' in row['comp_name'] else 'E'
        out_str += ' '*(6-len(anomoly_name)) + str(anomoly_name)
        out_str += ' '*(9-len(date)) + date
        out_str += ' '*(8-len("%.2f"%float(row['inter_lat']))) + "%.2f"%float(row['inter_lat'])
        out_str += ' '*(8-len("%.2f"%float(row['inter_lon']))) + "%.2f"%float(row['inter_lon'])
        out_str += ' '*(8-len("%.2f"%(float(row['strike'])))) + "%.2f"%(float(row['strike']))
        out_str += ' '*(7-len(phase_shift)) + phase_shift
        out_str += ' '*(11-len('10.000')) + '10.000' + '\n'
    if outfile==None: outfile = "maxtab.%s"%anomoly_name
    print("saving to %s"%outfile)
    out_file = open(outfile,'w+')
    out_file.write(out_str)
    out_file.close()

def create_max_file(deskew_df,srf,asf,outfile="deskew.max",aero_1s=10.,ship_1s=10.):
    #Open Out File
    with open(outfile,'w+') as fout:
        #Write Max Header
        fout.write("Max File\n")
        num_profiles = len(deskew_df[deskew_df["track_type"]=="aero"])/2 + len(deskew_df[deskew_df["track_type"]=="ship"])
        fout.write("80.0,0.0,1.0,0,0,%d,0,0,0,1,1.0\n"%(num_profiles))
        #Write Each Track
        deskew_df = calc_aei(deskew_df,srf,asf)
        for i,row in deskew_df.iterrows():
            if row["track_type"]=="ship":
                lat = row["inter_lat"]
                lon = row["inter_lon"]
                strike = row["strike"]
                aei = row["aei"]
                s1aei = ship_1s
            else:
                s1aei = aero_1s
                if "Ed" in row["comp_name"]:
                    other_comp = row["comp_name"].replace(".Ed.lp",".Vd.lp")
                    oth_row = deskew_df[deskew_df["comp_name"]==other_comp]
                    lat = (row["inter_lat"] + oth_row["inter_lat"])/2
                    lon = (row["inter_lon"] + oth_row["inter_lon"])/2
                    strike = (row["strike"] + oth_row["strike"])/2
                    aei = (row["aei"] + oth_row["aei"])/2
                elif "Hd" in row["comp_name"]:
                    other_comp = row["comp_name"].replace(".Hd.lp",".Vd.lp")
                    oth_row = deskew_df[deskew_df["comp_name"]==other_comp]
                    lat = (row["inter_lat"] + oth_row["inter_lat"])/2
                    lon = (row["inter_lon"] + oth_row["inter_lon"])/2
                    strike = (row["strike"] + oth_row["strike"])/2
                    aei = (row["aei"] + oth_row["aei"])/2
                elif "Vd" in row["comp_name"]: continue
                else: raise ValueError("Unknown aeromagnetic component type for %s"%str(row["comp_name"]))

            fout.write(row["comp_name"]+"\n")
            fout.write("%.2f,%.2f,%.2f,%.2f,%.2f\n"%(aei,s1aei,lat,lon,strike))
        #Write Fake Remanent Amp Factor to prevent singularity in old Max
        fout.write("Fake Amplitude\n")
        fout.write("1.0,0.1,0.0,180.0,90.0")

def generate_spreading_rate_model(spreading_rate_path):
    sr_file = open(spreading_rate_path,'r')
    sr_lines = sr_file.readlines()
    sr_dict,header = {},""
    for line in sr_lines:
        entries = line.split("\t")
        try:
            sr_dict[header].append(list(map(float,entries)))
        except (ValueError,KeyError) as e:
            header = line.strip("\n")
            if header=="": continue
            sr_dict[header] = []

    def spreading_rate_model(sz,age):
        if sz not in sr_dict.keys() and 'Default' in sr_dict.keys(): print("%s not in spreading rate model switching to default"%sz); sz = 'Default'
        for ub_age,sr in sr_dict[sz]:
            if float(age)<=ub_age: return sr
        print("age %.3f could not be found in the spreading rate model for spreading zone %s"%(float(age),sz))

    return spreading_rate_model,sr_dict.keys()

def generate_anomalous_skewness_model(anomalous_skewness_path):
    as_file = open(anomalous_skewness_path,'r')
    as_lists = [list(map(float,line.rstrip('\n').split('\t'))) for line in as_file.readlines()]

    def anomalous_skewness_model(sr):
        for i,(ub_sr,as_val) in enumerate(as_lists):
            if float(sr)<=ub_sr:
                if i==0: return as_val
                m = (as_val-as_lists[i-1][1])/(ub_sr-as_lists[i-1][0])
                b = as_val - m*ub_sr
                return m*sr+b
        return as_val
#        print("spreading rate %.3f could not be found in the anomalous skewness model"%float(sr))

    return anomalous_skewness_model

def get_shipmag_decimal_year(row,deg_e=.01):
    """
    takes a row (pandas series) of a deskew file which is of type ship and returns the decimal year for the
    intersection point. returns none if not found.
    """
    if row['track_type']!='ship':
        raise ValueError("get_shipmag_decimal_year can only run on shipmag data recieved data of type %s instead"%str(row['track_type']))
    data_file_path = os.path.join(row["data_dir"],row["comp_name"])
    data_df = utl.open_mag_file(data_file_path)
#    data_df = pd.read_csv(data_file_path,names=["dist","decimal_year","mag","lat","lon"],delim_whitespace=True)
    decimal_year=None
    for j,datarow in data_df.iterrows(): #iterate to find the distance associated with the current lat lon
        if (float(datarow['lat'])>=float(row['inter_lat'])-deg_e and \
          float(datarow['lat'])<=float(row['inter_lat'])+deg_e) and \
          (utl.convert_to_0_360(datarow['lon'])>=utl.convert_to_0_360(row['inter_lon'])-deg_e and \
          utl.convert_to_0_360(datarow['lon'])<=utl.convert_to_0_360(row['inter_lon'])+deg_e):
            decimal_year=float(datarow['dec_year']); break
    return decimal_year

def new_get_shipmag_decimal_year(row,deg_e=.01):
    """
    takes a row (pandas series) of a deskew file which is of type ship and returns the decimal year for the
    intersection point. returns none if not found.
    """
    if row['track_type']!='ship':
        raise ValueError("get_shipmag_decimal_year can only run on shipmag data recieved data of type %s instead"%str(row['track_type']))
    data_file_path = os.path.join(row["data_dir"],row["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","decimal_year","mag","lat","lon"],delim_whitespace=True)
    decimal_year,sum_sq=None,1e9
    for j,datarow in data_df.iterrows(): #iterate to find the distance associated with the current lat lon
        next_sum_sq = (float(row['inter_lat'])-float(datarow['lat']))**2 + (utl.convert_to_0_360(row['inter_lon'])-utl.convert_to_0_360(datarow['lon']))**2
        if sum_sq<next_sum_sq: break
        elif (float(datarow['lat'])>=float(row['inter_lat'])-deg_e and \
          float(datarow['lat'])<=float(row['inter_lat'])+deg_e) and \
          (utl.convert_to_0_360(datarow['lon'])>=utl.convert_to_0_360(row['inter_lon'])-deg_e and \
          utl.convert_to_0_360(datarow['lon'])<=utl.convert_to_0_360(row['inter_lon'])+deg_e):
            sum_sq = next_sum_sq
            decimal_year=float(datarow['dec_year'])
    return decimal_year

def reduce_to_pole(deskew_path, pole_lon, pole_lat, spreading_rate_path=None, anomalous_skewness_model_path=None):

    asf,srf,sz_list = get_asf_srf(spreading_rate_path,anomalous_skewness_model_path)

    deskew_df = filter_deskew_and_calc_aei(deskew_path)

    print("reducing to pole - lat: %.3f, lon: %.3f"%(pole_lat,pole_lon))

    deskew_df = reduce_dsk_df_to_pole(deskew_df, pole_lon, pole_lat, asf, srf)

    old_results_dir = deskew_df['results_dir'].iloc[0]
    new_results_dir = os.path.join(old_results_dir,"pole_%.0f_%.0f_results"%(pole_lon,pole_lat))
    utl.check_dir(new_results_dir)
    deskew_df['results_dir'] = new_results_dir

#    reduced_deskew_df = deskew_df[['comp_name','phase_shift','step','rel_amp','age_min','age_max','inter_lat','inter_lon','strike','data_dir','results_dir','track_type','sz_name','r','g','b']]

    out_path = os.path.join(os.path.dirname(deskew_path),"pole_%.0f_%.0f.deskew"%(pole_lon,pole_lat))
    print("writing to %s"%out_path)
    utl.write_deskew_file(out_path,deskew_df)

def reduce_dsk_df_to_pole(dsk_df, pole_lon, pole_lat, asf, srf):

    deskew_df = dsk_df.copy()

    if "phase_shift" not in deskew_df.columns: deskew_df["phase_shift"] = 0.
    if "rel_amp" not in deskew_df.columns: deskew_df["rel_amp"] = 0.

    for i,row in deskew_df.sort_values("comp_name",ascending=False).iterrows():
        e_r,reduced_skewness,rel_reduced_amplitude = reduce_dsk_row_to_pole(row, pole_lon, pole_lat, asf, srf)

        deskew_df.at[i,'phase_shift'] = round(reduced_skewness,3)
        deskew_df.at[i,'aei'] = round(e_r,3)
#        if "Ed.lp" in row["comp_name"]:
#            rel_reduced_amplitude = deskew_df[deskew_df["comp_name"]==row["comp_name"].replace("Ed.lp","Vd.lp")].iloc[0]["rel_amp"]
#        elif "Hd.lp" in row["comp_name"]:
#            rel_reduced_amplitude = deskew_df[deskew_df["comp_name"]==row["comp_name"].replace("Hd.lp","Vd.lp")].iloc[0]["rel_amp"]
#        deskew_df.at[i,'rel_amp'] = round(rel_reduced_amplitude,3)

    return deskew_df

def reduce_dsk_row_to_pole(row, pole_lon, pole_lat, asf, srf):
        #taken from Lin's Matlab Code (doesn't seem to work)
#        geodes_inv_dic = Geodesic.WGS84.Inverse(float(row['inter_lat']),float(row['inter_lon']),pole_lat,pole_lon)
#        strike,arc_len_pole,az_pole = float(row['strike']),geodes_inv_dic['a12'],geodes_inv_dic['azi1']
#        I = np.arctan(2*np.tan(np.deg2rad(90-arc_len_pole)))
#        a = az_pole+180-strike
#        e = utl.wrap_180_180(np.degrees(np.arctan2(np.tan(np.deg2rad(I)),np.sin(np.deg2rad(a)))))
#        reduced_skewness = utl.wrap_180_180(float(row['phase_shift']) - (e - float(row['ei'])))

        rad_inter_lat,rad_inter_lon = np.deg2rad(float(row['inter_lat'])),np.deg2rad(float(row['inter_lon']))
        rad_pole_lat,rad_pole_lon = np.deg2rad(pole_lat),np.deg2rad(pole_lon)
        rad_strike = np.deg2rad(float(row['strike']))

        g = np.cos(rad_pole_lat)*np.cos(rad_inter_lat)*np.cos(rad_pole_lon-rad_inter_lon) + np.sin(rad_pole_lat)*np.sin(rad_inter_lat) #Angle between pole and site

        I = np.arctan2((2*g),np.sqrt(1-g**2))

        Dnum = np.cos(rad_pole_lat)*np.sin(rad_inter_lon-rad_pole_lon)
        Dden = np.cos(rad_inter_lat)*np.sin(rad_pole_lat) - np.sin(rad_inter_lat)*np.cos(rad_pole_lat)*np.cos(rad_inter_lon-rad_pole_lon)
        D = -np.arctan2(Dnum,Dden)

        e_r = np.degrees(np.arctan2(np.tan(I),np.sin(rad_strike-D)))

        anom_skew = asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2))

        reduced_skewness = utl.wrap_0_360(180. - float(row['ei']) - e_r + anom_skew)
        rel_reduced_amplitude = (np.sin(np.deg2rad(float(e_r)))/np.sin(np.deg2rad(row['ei'])))*np.sqrt(1+3*(np.cos(np.deg2rad(90-row["inter_lat"])))**2)

        return e_r,reduced_skewness,rel_reduced_amplitude

def create_deskew_file(chron_name,results_directory,age_min,age_max,data_directory='.',phase_shift=180,step=60):
    cut_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
    cut_tracks_file = open(cut_tracks_path,'r')
    cut_tracks = cut_tracks_file.readlines()
    cut_tracks_file.close()
    track_sz_and_inters = [track.split('\t') for track in cut_tracks]

    track_sz_and_inters.sort()

#    dout = {'comp_name':[],'phase_shift':[],'step':[],'inter_lat':[],'inter_lon':[],'strike':[],'age_min':[],'age_max':[],'track_type':[],'data_dir':[],'results_dir':[]}
#    dfout = pd.DataFrame(columns=['comp_name','phase_shift','step','inter_lat','inter_lon','strike','age_min','age_max','track_type','data_dir','results_dir'])
    out_str="comp_name\tphase_shift\tstep\tage_min\tage_max\tinter_lat\tinter_lon\tstrike\tdata_dir\tresults_dir\ttrack_type\tsz_name\tr\tg\tb\n"
    colors,i,j = [(0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89), (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236), (255, 188, 121), (207, 207, 207)],0,1
    for i in range(len(colors)):
        r, g, b = colors[i]  
        colors[i] = (r / 255., g / 255., b / 255.)
    sz_to_color,sz_to_name = {},{}
    for track,sz,inter in track_sz_and_inters:

        idx = utl.read_idx_from_string(inter)
        if 'aero' in track:
            Ef = os.path.basename(track) + '.Ed.lp'
            Vf = os.path.basename(track) + '.Vd.lp'
            comps = [Ef,Vf]
            track_type = 'aero'
        elif 'ship' in track:
            comps = [os.path.basename(track) + '.lp']
            track_type = 'ship'

        data_dir = os.path.split(track)[0]

        azsz_path = track[:-3]+'_'+track[-2:]+'_'+os.path.basename(sz)[:-4]+'.azszs'
        try:
            azsz_df = pd.read_csv(azsz_path,sep=' ')
            strike = azsz_df['strike'][0]
        except:
            strike = 180

        if sz not in sz_to_color:
            if i >= len(colors)-1: print("more spreading zones than colors, looping"); i = 0
            sz_to_color[sz] = colors[i]
            sz_to_name[sz] = "Spreading Zone %d"%j
            i += 1; j += 1

        for f in comps:
            out_str += f + "\t"
            out_str += str(phase_shift) + "\t"
            out_str += str(step) + "\t"
            out_str += str(age_min) + "\t"
            out_str += str(age_max) + "\t"
            out_str += str(idx[0][0][1]) + "\t"
            out_str += str(idx[0][0][0]) + "\t"
            out_str += str(strike) + "\t"
            out_str += str(data_dir) + "\t"
            out_str += str(results_directory) + "\t"
            out_str += track_type + "\t"
            out_str += str(sz_to_name[sz]) + "\t"
            out_str += str(sz_to_color[sz][0]) + "\t"
            out_str += str(sz_to_color[sz][1]) + "\t"
            out_str += str(sz_to_color[sz][2]) + "\n"

    utl.check_dir(data_directory)
    fout = open(os.path.join(data_directory,chron_name+'.deskew'), "w+")
    fout.write(out_str)
    fout.close()

def create_matlab_datatable(deskew_path,deg_e=.0001):
    """Creates Datatable"""
    deskew_df = filter_deskew_and_calc_aei(deskew_path)

    deskew_df['r'],deskew_df['g'],deskew_df['b'] = 0,0,0 #overwriting color because this is only used for the deskew plots and colors are annoying there

    dt = deskew_df[["comp_name","data_dir","inter_lat","inter_lon","phase_shift","aei","strike","r","g","b"]]

    dt_path = "%s.dt.csv"%os.path.basename(deskew_path).split('.')[0]
    dt.to_csv(dt_path,sep="\t",index=False,header=False)

    return dt_path,deskew_df["age_min"].iloc[0],deskew_df["age_max"].iloc[0],deskew_df["results_dir"].iloc[0]

def create_spreading_rate_file(spreading_rate_picks_path, ages_path):
    spreading_rate_picks = pd.read_csv(spreading_rate_picks_path,sep='\t',header=0,index_col=0)
    ages_df = pd.read_csv(ages_path,sep='\t',header=0,index_col=0)

    first_anom = ages_df.index[len(ages_df)-1]
    
    sr_dict = {}
    for track in spreading_rate_picks.columns:
        sr_dict[track] = {}
        track_picks = spreading_rate_picks[spreading_rate_picks[track].notnull()][track]
        prev_anom = track_picks.index[0]
        
        # Find first anomaly in this spreading zone
        if (ages_df.index.get_loc(prev_anom) < ages_df.index.get_loc(first_anom)): first_anom = prev_anom
        
        for anom,pick in track_picks.iloc[1:].iteritems():
            width = track_picks[prev_anom] - pick
            duration = ages_df['base'][anom] - ages_df['base'][prev_anom]
            sr_dict[track][anom] = width/duration
            prev_anom = anom
    sr_df = pd.DataFrame(sr_dict, index=spreading_rate_picks.index[1:])
    
    # Save file with spreading rate estimates for individual tracks
    #import pdb; pdb.set_trace()
    if not os.path.isdir(os.path.join(os.path.dirname(spreading_rate_picks_path),'profile_sr')): os.makedirs(os.path.join(os.path.dirname(spreading_rate_picks_path),'profile_sr'))
    track_path = os.path.join(os.path.dirname(spreading_rate_picks_path),'profile_sr',os.path.basename(spreading_rate_picks_path).split('.')[0] + '_profiles.sr')
    sr_df.to_csv(track_path,sep='\t')

    stats_dict = {'mean':{},'median':{},'std':{},'n':{},'age_min':{},'age_max':{}}
    # The commented line below always starts the age domain at the youngest anomaly. This may not be accurate if picks don't go that far (e.g. CL-CL).
    #prev_anom = ages_df.index[list(ages_df.index).index(sr_df.index[0])-1]
    prev_anom = first_anom
    for anom,picks in sr_df.iterrows():
        pick_list = picks[picks.notnull()].tolist()
        if not pick_list: print("no data to average for %s, skipping"%anom); continue
        stats_dict['mean'][anom] = "%.3f"%stat.mean(pick_list)
        stats_dict['median'][anom] = "%.3f"%stat.median(pick_list)
        if len(pick_list)>1:
            stats_dict['std'][anom] = "%.3f"%stat.stdev(pick_list)
        else:
            stats_dict['std'][anom] = "%.3f"%0.000
        stats_dict['n'][anom] = "%.3f"%len(pick_list)
        stats_dict['age_min'][anom] = "%.3f"%ages_df['base'][prev_anom]
        stats_dict['age_max'][anom] = "%.3f"%ages_df['base'][anom]
        prev_anom = anom
    stats_df = pd.DataFrame(stats_dict, index=sr_df.index)
    stats_df = stats_df[stats_df['mean'].notnull()]

    out_path = os.path.join(os.path.dirname(spreading_rate_picks_path),os.path.basename(spreading_rate_picks_path).split('.')[0] + '.sr')
    stats_df.to_csv(out_path,sep='\t')  

    return stats_df

def get_lon_lat_from_plot_picks_and_deskew_file(deskew_path,spreading_rate_picks_path):
    deskew_df = utl.open_deskew_file(deskew_path)
    spreading_rate_picks = pd.read_csv(spreading_rate_picks_path,sep='\t',header=0,index_col=0)

    iso_dict = {}
    lats,lons = [],[]
    for track in spreading_rate_picks.columns:
        iso_dict[track] = {}
        track_picks = spreading_rate_picks[track]
        track_picks_tups = tuple(track_picks.groupby((track_picks.isnull()!=track_picks.isnull().shift()).cumsum()))
        for i,track_picks in track_picks_tups:
            if track_picks[track_picks.notnull()].empty: continue
            prev_anom = track_picks.index[0]
            for anom,pick in track_picks.iloc[1:].iteritems():
                iso_dict[track][anom] = {}
                iso_dist = (pick+track_picks[prev_anom])/2
                drow = deskew_df[deskew_df['comp_name'].str.rstrip('.Ed .Vd .lp') == track]
                if drow.empty: print('problem getting deskew data from spreading rate picks for %s, check track names for now skipping'%track); return pd.DataFrame(),0,0
                drow = drow.iloc[0] #make sure it is the first value and that it is a series
                iso_lon,iso_lat,iso_dist_on_track = get_lon_lat_from_plot_pick(drow,iso_dist)
                iso_dict[track][anom]['lon'] = iso_lon
                iso_dict[track][anom]['lat'] = iso_lat
                lons.append(iso_lon)
                lats.append(iso_lat)
                prev_anom = anom
    iso_df = pd.DataFrame(iso_dict)
    average_lat = sum(lats)/len(lats)
    average_lon = sum(lons)/len(lons)

    # Initialize empty dataframe to hold isochron picks
    pick_df = pd.DataFrame(index=iso_df.columns, columns=['lon','lat'])
    for anom in iso_df.index:
        for track in iso_df.columns:
            pick_df.loc[track] = iso_df[track][anom]
        # Flip columns so latitude comes before longitude, because we aren't savages, Kevin
        pick_df = pick_df.reindex(columns=['lat','lon'])
        # Drop rows without data
        pick_df = pick_df.dropna()
        # Define path and filename
        picks_path = os.path.join(os.path.dirname(spreading_rate_picks_path),os.path.basename(spreading_rate_picks_path).split('.')[0] + '_' + anom + '.txt')
        #import pdb; pdb.set_trace()
        pick_df.to_csv(picks_path,sep='\t')

    return iso_df,average_lon,average_lat

def make_spreading_rate_model(sr_paths):
    out_string = ""
    for sr_path in sr_paths:
        try: sr_df = pd.read_csv(sr_path,sep='\t')
        except IOError: raise IOError("There can be no other flags or values besides .sr file paths after the -srm flag, was given %s as a .sr file path and it either doesn't exist or is completely bananas"%sr_path)
        out_string += '------Please rename with the sz for ' + sr_path + '------\n'
        for i,row in sr_df.iterrows():
            age_max = "%.3f"%float(row['age_max'])
            if i == sr_df.index[-1]: age_max = "1e10"
            out_string += age_max + '\t' + "%.3f"%float(row['mean']) + '\n'
    out_file = open('new_spreading_rate_model.txt','w+')
    out_file.write(out_string)
    out_file.close()

def save_iso_picks(deskew_path,srp_paths):
    anom_out_str = {}
    for srp_path in srp_paths:
        iso_df,avg_lon,avg_lat = get_lon_lat_from_plot_picks_and_deskew_file(deskew_path,srp_path)
        for anom,row in iso_df.iterrows():
            if anom in anom_out_str.keys(): anom_out_str[anom] += "> %s\n"%srp_path
            else: anom_out_str[anom] = "> %s\n"%srp_path
            row = row[row.notnull()]
            for lon_lat_dict in sorted(row.tolist(),key=lambda x,y={'lat':0}: x['lat']-y['lat']):
                anom_out_str[anom] += "%.3f\t%.3f\n"%(float(lon_lat_dict['lon']),float(lon_lat_dict['lat']))
    for key in anom_out_str.keys():
        out_dir = "new_isochron_picks"
        utl.check_dir(out_dir)
        out_file = "%s_chron.txt"%key
        outf = os.path.join(out_dir,out_file)
        f_out = open(outf,"+w")
        f_out.write(anom_out_str[key])
        f_out.close()

def read_and_fit_ridge_data(ridge_loc_path):
    ridge_loc_file = open(ridge_loc_path,'r')
    lines = ridge_loc_file.readlines()

    ridge_loc_dict,current_head = {},''
    for line in lines:
        if line.startswith('>'):
            current_head = line.strip('> \n')
            ridge_loc_dict[current_head] = []
        else:
            ridge_loc_dict[current_head].append(list(map(float,line.strip('\n').split('\t'))))

    for sz in list(ridge_loc_dict.keys()):
        if ridge_loc_dict[sz]==[]: ridge_loc_dict.pop(sz); continue
        ar = np.array(ridge_loc_dict[sz])
        dists = [Geodesic.WGS84.Inverse(ar[0,1],utl.convert_to_0_360(ar[0,0]),lat,utl.convert_to_0_360(lon))['s12']/1000 for lat,lon in zip(ar[:,1],utl.convert_to_0_360(ar[:,0]))]
        new_dists = np.arange(0,dists[-1],10)
        ridge_lons = np.interp(new_dists,dists,utl.convert_to_0_360(ar[:,0]))
        ridge_min_lon = min(ridge_lons)
        ridge_max_lon = max(ridge_lons)
        ridge_lats = np.interp(new_dists,dists,ar[:,1])
        ridge_min_lat = min(ridge_lats)
        ridge_max_lat = max(ridge_lats)
        ridge_lon_lats = [[lon,lat] for lon,lat in zip(ridge_lons,ridge_lats)]
        ridge_loc_dict[sz] = ridge_lon_lats

    def get_ridge_loc(sz,track_lon_lats):
        if sz not in ridge_loc_dict.keys(): print("sz not found when looking for ridge location, was given spreading zone %s but only had options %s"%(str(sz),str(ridge_loc_dict.keys()))); return None,None
        idx = gp.intersect_bf(track_lon_lats, ridge_loc_dict[sz],e=.5)
        if idx == [None,None]: print("could not calculate intersept"); return None,None
        else: return idx[0][0]

    return get_ridge_loc

def read_and_fit_fz_data(fz_directory=os.path.join('..','raw_data','fracture_zones')):
    fzs = glob.glob(os.path.join(fz_directory,'*'))
    lfz = []
    for fz in fzs:
        fzdf = pd.read_csv(fz,sep='\t')
        dists = [Geodesic.WGS84.Inverse(fzdf['Latitude'][0],utl.convert_to_0_360(fzdf['Longitude'][0]),lat,utl.convert_to_0_360(lon))['s12']/1000 for lat,lon in zip(fzdf['Latitude'],fzdf['Longitude'])]
        new_dists = np.arange(0,dists[-1],10)
        fz_lons = np.interp(new_dists,dists,utl.convert_to_0_360(fzdf['Longitude']))
        fz_lats = np.interp(new_dists,dists,fzdf['Latitude'])
        fz_lon_lats = [[lon,lat] for lon,lat in zip(fz_lons,fz_lats)]
        lfz.append(fz_lon_lats)

    def get_fz_loc(track_lon_lats,e=.5):
        inters = []
        for fz in lfz:
            idx = intersect_bf(track_lon_lats, fz,e=e)
            if idx != [None,None]: inters.append(idx[0][0])
        return inters

    return get_fz_loc

def find_fz_crossings(deskew_path,fz_directory=os.path.join('..','raw_data','fracture_zones')):
    deskew_df = utl.open_deskew_file(deskew_path)
    get_fz_loc = read_and_fit_fz_data(fz_directory)

    fz_inter_dict = {}
    for i,row in deskew_df.iterrows():
        data_path = os.path.join(row['data_dir'],row['comp_name'])
        data_df = utl.open_mag_file(data_path)
        track_lon_lats = [[utl.convert_to_0_360(lon),lat] for lon,lat in zip(data_df['lon'],data_df['lat'])]
        inters = get_fz_loc(track_lon_lats)
        if inters != []: fz_inter_dict[row['comp_name']] = inters

    fz_inter_df = pd.DataFrame({'inters':fz_inter_dict})
    fz_inter_df.to_csv('fz_intercepts.txt',sep='\t')

def update_useable_tracks_from_deskew(deskew_path,useable_track_path):
    useable_df = pd.read_csv(useable_track_path, sep='\t', header=None)
    useable_df['tracks'] = list(map(os.path.basename, useable_df[0].tolist()))
    deskew_df = utl.open_deskew_file(deskew_path)
    useable_tracks = list(map(lambda x: x.rstrip('.Ed .lp .Vd'), deskew_df['comp_name'].tolist()))
    new_useable_df = useable_df[useable_df['tracks'].isin(useable_tracks)][[0,1,2]]
    directory = os.path.dirname(useable_track_path)
    new_useable_track_filename = 'new_' + os.path.basename(useable_track_path)
    out_path = os.path.join(directory,new_useable_track_filename)
    new_useable_df.to_csv(out_path, sep='\t', index=False, header=False)

def create_deskewed_data_file(deskew_path):
    #read deskew file
    deskew_df = utl.open_deskew_file(deskew_path)

    #iterate mag files
    for i,row in deskew_df.iterrows():
        #read mag files
        data_path = os.path.join(row['data_dir'],row['comp_name'])
        data_df = utl.open_mag_file(data_path)
        #deskew mag data
        data_df['deskewed_mag'] = phase_shift_data(data_df['mag'],float(row['phase_shift']))
        #save deskewed mag data as $DATAFILE.deskew
        print("writing %s"%(data_path+'.deskewed'))
        data_df[['lon','lat','deskewed_mag']].to_csv(data_path+'.deskewed',sep=',',header=False,index=False)

def config_synthetic(synth_config_path,timescale_path,spreading_rate_path,length=4096,buff=.2):
    sconf = pd.read_csv(synth_config_path,sep='\t',header=0)
    synth_and_domains = []
    for i,row in sconf.iterrows():
        synth_and_domains.append(make_synthetic(*row.values,timescale_path,spreading_rate_path,length=length,buff=buff))
    return synth_and_domains

def make_synthetic(age_min,age_max,layer_depth,layer_thickness,layer_mag,azi,rd,ri,ad,ai,fix_sta,fix_end,twf,timescale_path,spreading_rate=None,sz_name=None,spreading_rate_path=None,length=4096):

    tdf = pd.read_csv(timescale_path,sep='\t',header=0)
    if spreading_rate_path!=None and sz_name!=None:
        srf = generate_spreading_rate_model(spreading_rate_path)[0]
    elif spreading_rate!=None:
        srf = lambda x,y: spreading_rate
    else:
        raise ValueError("Either a constant spreading rate must be given or a spreading rate model file and a spreading zone name in that file")

    #Construct Domains and Square Wave
    age_step = (age_max-age_min)/(length-1) #time sampling rate
    #Start total distance at min_distance from central
    if age_min>0: min_dis = sum(map(lambda x: age_step*srf(sz_name,x),np.arange(0,age_min,age_step)))
    else: min_dis = -sum(map(lambda x: age_step*srf(sz_name,x),np.arange(age_min,0,age_step)))
    tot_dis,mag_sig,time_domain,dis_domain = 0,[],np.linspace(age_min,age_max,length),[]
    for age in time_domain:
        tot_dis += age_step*srf(sz_name,age)
        dis_domain.append(tot_dis+min_dis)
        try: idx = tdf[(abs(round(age,2))>=tdf["top"]) & (abs(round(age,2))<=tdf["base"])].index[0]
        except: raise IndexError("Age bounds of %.3f, %.3f out of bounds of age model %s"%(age_min,age_max,timescale))
        if idx%2==0: pol = layer_mag
        else: pol = -layer_mag
        mag_sig.append(pol) #Square Wave
    samp_dis = (tot_dis)/(length-1) #Distance sampling rate

    #Pad Signal before processing
    # ibegin = (NDIM-1)*(.15/1.3)+1
    # iend = int((NDIM-1)*(1.15/1.3) + 0.5)
    buff = int((2**(int(np.log2(length))+1)-length)/2)
    if not fix_sta: lbuff = mag_sig[0]*np.ones([int(buff)])
    else: lbuff = np.zeros([int(buff)])
    if not fix_sta: rbuff = mag_sig[-1]*np.ones([int(buff)])
    else: rbuff = np.zeros([int(buff)])

    mag_sig = np.hstack((lbuff,mag_sig,rbuff))

    #Calc effective inclination
    a_ei = np.rad2deg(np.arctan2(np.tan(np.deg2rad(ai)),np.cos(np.deg2rad(ad-azi)))) #Apparent Effective Inc
    r_ei = np.rad2deg(np.arctan2(np.tan(np.deg2rad(ri)),np.cos(np.deg2rad(rd-azi)))) #Remanent Effective Inc
    rtheta = np.deg2rad(a_ei+r_ei+180) #Remanent Phase shift
    afac = (np.sin(np.deg2rad(ai))*np.sin(np.deg2rad(ri)))/(np.sin(np.deg2rad(a_ei))*np.sin(np.deg2rad(r_ei))) #Remanent Amplitude Factor

    ft_mag_sig = np.fft.fft(mag_sig) #Transform to freq domain

    dis_nyquist = np.pi/samp_dis #We're really concerned with distance domain so consider sampling in that space
    #This is the important step where we convolve with the earth's response
    fte_mag_sig = mag_earth_filter(ft_mag_sig,length+2*buff,layer_depth,layer_depth+layer_thickness,rtheta,afac,dis_nyquist)

    if twf: #Guess what this does.
        fte_mag_sig = transition_width_filter(fte_mag_sig,length+2*buff,dis_nyquist,twf/4)

    synth = np.fft.ifft(fte_mag_sig)

    return synth[int(buff):-int(buff)],dis_domain,samp_dis
#    return synth,dis_domain,samp_dis

def mag_earth_filter(ft_mag_sig,length,top,bot,rtheta,afac,dis_nyquist):
    z = np.cos(rtheta)+1j*np.sin(rtheta)
    ft_mag_sig[0] = np.imag(ft_mag_sig[0])
    sinc = dis_nyquist/(length/2)
    for i,j in zip(range(int(length/2)),range(int(length-1),int(length/2),-1)):
        sdis = i*sinc
        efac = 2*np.pi*(np.exp(-top*sdis)-np.exp(-bot*sdis))*afac
        ft_mag_sig[i] = efac*z*ft_mag_sig[i]
        ft_mag_sig[j] = efac*np.conjugate(z)*ft_mag_sig[j]
    return ft_mag_sig


def transition_width_filter(sig,length,nyquest,sigma):
    for i,j in zip(range(int(length/2)),range(int(length-1),int(length/2),-1)):
        sdis = i*(nyquest/(length/2))
        gfac = np.exp(-(np.pi*sigma*sdis)**2)
        sig[i] = sig[i]*gfac
        sig[j] = sig[j]*gfac
    return sig


def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation. 
    Shifted data filled with NaNs 
    
    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects or objects castable to pd. Series of equal length
    Returns
    ----------
    crosscorr : float
    """
    datax,datay = pd.Series(datax),pd.Series(datay)
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else: 
        return datax.corr(datay.shift(lag))

def butter_lowpass(highcut, fs, order=5):
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = butter(order, high, btype='lowpass')
    return b, a


def butter_lowpass_filter(data, highcut, fs, order=5):
    b, a = butter_lowpass(highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def auto_dsk(dsk_row,synth,bounds,conv_limit=0,conv_bounds=[None,None],phase_args=(0.,360.,1.),highcut=0.,order=3):
    """
    Returns the maximum likelihood phase shift to deskew the data to match a provided synthetic given a bounds
    on a window to match.

    Parameters
    ----------

    dsk_row : Pandas.Series
        Single row of a deskew file with valid path to data file
    synth : list
        This should be the output of the make synthetic function it needs to contain three elements
            0) an array of the synthetic magnetic anomalies
            1) an array of the distance coordinates of the points in 0, should be of equal length to 0, MUST be in
               the same coordinate system as the profile provided in dsk_row!!! Which it may not by default.
            2) the distance resolution of the synthetic in 0 and 1
    bounds : list of floats
        Has two elements which corespond to the left and right bounds of the window
    conv_limit : float, optional
        Weather or not to realign the anomaly each phase shift using a time lagged convolution method which
        increases runtime significantly but can also increase accuracy. This argument should be a positve
        float which corresponds to the amount of +- shift the anomaly is allowed to move used otherwise it should
        be 0 to not use the shift method (Default: 0, which implies not to use method).
    conv_bounds : list of 2 floats, optional
        The left and right boundary in the distance domain to use to time lag convolve the synthetic and the filtered
        data signal. Thus 300 km of signal can be convolved but only the 10 km of motion allowed to pin down the 
        crossing location. (Default: [None,None], which implies conv_bounds=bounds)
    phase_args : tuple or other unpackable sequence, optional
        Arguments to np.arange which define the phases searched in the minimization. (Default: (0.,360.,1.) which
        implies a search of the entire parameter space of phases at 1 degree resolution)
    highcut : float, optional
        The upper cutoff frequency to filter the data by in order to remove any topographic anomalies in the data.
        This value should be between 0 and Nyquest of the synthetic which MUST be regularly sampled like those
        returned by make_synthetic. The data is up or down sampled to the synthetic before filtering. (Default:
        0 which implies not to filter the data)
    order : int, optional
        The order of the lowpass butterworth filter to apply to the data.

    Returns
    ----------

    best_phase : float
        The maximum liklehood phase shift to match the data to the synthetic
    best_shift : float
        the maximum likelihood shift for the best_phase which aligned the two anomalies
    phase_func : Numpy.NdArray
        The summed phase asynchrony between the data and the synthetic as a function of phase shift (best_phase is
        the global minimum of this function)
    best_shifts : Numpy.NdArray
        the maximum likelihood shift as a function of the phase shift
    """

    #Unpack Arguments
    dage = dsk_row["age_max"]-dsk_row["age_min"]
    phases = np.arange(*phase_args)
    left_bound,right_bound = bounds
    synth_mag = np.array(synth[0])
    synth_dis = np.array(synth[1])
    ddis = synth[2]

    data_path = os.path.join(dsk_row["data_dir"],dsk_row["comp_name"])
    data_df = utl.open_mag_file(data_path)
    projected_distances = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),(180+dsk_row['strike'])%360)
    if conv_limit: #create the fully interpolated profile for convolution

        #create the shortened synthetic for the time lagged convolution
        if isinstance(conv_bounds[0],type(None)): conv_bounds[0] = bounds[0]
        if isinstance(conv_bounds[1],type(None)): conv_bounds[1] = bounds[1]
        left_idx = np.argmin(np.abs(synth_dis - conv_bounds[0]))
        right_idx = np.argmin(np.abs(synth_dis - conv_bounds[1]))
        right_idx,left_idx = max([right_idx,left_idx]),min([right_idx,left_idx])
        conv_synth,conv_synth_dis = synth_mag[left_idx:right_idx],synth_dis[left_idx:right_idx]

        if np.any(np.diff(projected_distances["dist"])<0): #redefine to the right because interp dumbs
            mag = data_df["mag"].to_numpy()[::-1]
            mag_dis = projected_distances["dist"].to_numpy()[::-1]
        full_imag = np.interp(conv_synth_dis,mag_dis,mag)
        if highcut: full_fimag = butter_lowpass_filter(full_imag,highcut=highcut,fs=1/ddis,order=order)
        else: full_fimag = full_imag

    #trim to only window of relivence
    left_idx = np.argmin(np.abs(synth_dis - left_bound))
    right_idx = np.argmin(np.abs(synth_dis - right_bound))
    right_idx,left_idx = max([right_idx,left_idx]),min([right_idx,left_idx])
    tsynth_mag = synth_mag[left_idx:right_idx]
    tsynth_dis = synth_dis[left_idx:right_idx]
    N = len(tsynth_mag) #because this is easier and regularly sampled plus the user can set it simply
    al2 = np.angle(hilbert(np.real(tsynth_mag),N),deg=False)

    best_shifts = [] #record best shifts as function of phase shift
    phase_async_func = [] #record summed phase asynchrony as a function of phase shift
    for i,phase in enumerate(phases):
        shifted_mag = phase_shift_data(data_df["mag"],phase)

        if conv_limit: #DON'T YOU KNOW WE'RE GONNAAAA DOOOOOOO THE COOONVOLUTIOOOON!!!
            shifted_full_fimag = phase_shift_data(full_fimag,phase)
            correlation_func = np.abs(np.convolve(shifted_full_fimag,conv_synth,"full"))
            correlation_func = correlation_func[int(len(conv_synth)-conv_limit/ddis+.5):int(len(conv_synth)+conv_limit/ddis+.5)]

            best_shift = ddis*(len(correlation_func)/2-np.argmax(correlation_func))/2

        else: best_shift = 0.

        #trim the data to the right segments
        left_idx = np.argmin(np.abs(projected_distances["dist"] - left_bound + best_shift))
        right_idx = np.argmin(np.abs(projected_distances["dist"]- right_bound + best_shift))
        right_idx,left_idx = max([right_idx,left_idx]),min([right_idx,left_idx])
        tproj_dist = projected_distances["dist"][left_idx:right_idx] + best_shift
        tshifted_mag = shifted_mag[left_idx:right_idx]

        #numpy.interp only works for monotonic increasing independent variable data
        if np.any(np.diff(tproj_dist)<0): itshifted_mag = np.interp(-tsynth_dis,-tproj_dist,tshifted_mag)
        else: itshifted_mag = np.interp(tsynth_dis,tproj_dist,tshifted_mag)
        if highcut: fitshifted_mag = butter_lowpass_filter(itshifted_mag,highcut=highcut,fs=1/ddis,order=order)
        else: fitshifted_mag = itshifted_mag

        al1 = np.angle(hilbert(fitshifted_mag,N),deg=False)
        phase_asynchrony = np.sin(np.abs(al1-al2)/2) #shouldn't go negative but...just in case
        best_shifts.append(best_shift)
        phase_async_func.append(phase_asynchrony.sum())

    best_idx = np.argmin(phase_async_func)

    return phases[best_idx],best_shifts[best_idx],phase_async_func,best_shifts


def find_best_phase_match(s0,s1,phase_args=(0.,360.,1.)):

    phases = np.arange(*phase_args)
    N = len(s1) #because this is easier and regularly sampled plus the user can set it simply
    al1 = np.angle(hilbert(np.real(s1),N),deg=False)

    best_shifts = [] #record best shifts as function of phase shift
    phase_async_func = [] #record summed phase asynchrony as a function of phase shift
    for i,phase in enumerate(phases):
        shifted_mag = phase_shift_data(s0,phase)

        al0 = np.angle(hilbert(shifted_mag,N),deg=False)
        phase_asynchrony = np.sin(np.abs(al0-al1)/2) #shouldn't go negative but...just in case
        best_shifts.append(0.)
        phase_async_func.append(phase_asynchrony.sum())

    best_idx = np.argmin(phase_async_func)

    return phases[best_idx],best_shifts[best_idx],phase_async_func,best_shifts
