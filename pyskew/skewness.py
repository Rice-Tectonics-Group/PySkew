import os, glob
import shutil
import subprocess
import pandas as pd
import numpy as np
import statistics as stat
import pmagpy.ipmag as ipmag
from multiprocessing import Process
from geographiclib.geodesic import Geodesic
from .geographic_preprocessing import *
from .utilities import *

def filter_deskew_and_calc_aei(deskew_path,spreading_rate_model_path=None,anomalous_skewness_model_path=None):
    """Creates Datatable"""

    asf,srf,sz_list = get_asf_srf(spreading_rate_model_path,anomalous_skewness_model_path)

    deskew_df = pd.read_csv(deskew_path,sep='\t',dtype=str)

    deskew_df = deskew_df[~deskew_df["comp_name"].str.startswith('#')]

    deskew_df["aei"] = [wrap_180_180(180-wrap_180_180(row['phase_shift'])-(90 if ".Vd." in row['comp_name'] else 0) + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2))) for i,row in deskew_df.iterrows()]

    deskew_df["ei"] = [(90 if ".Vd." in comp else 0) for ps,comp in zip(deskew_df["phase_shift"],deskew_df["comp_name"])]

    for i,row in deskew_df[deskew_df['track_type']=='ship'].iterrows():
        decimal_year = get_shipmag_decimal_year(row)
        if decimal_year==None: print("Intersection point could not be found in data file so IGRF could not be calculated and aei could not be found please check your data or increase degree error of create_matlab_datatable, skipping %s"%row['comp_name']); continue
        igrf = ipmag.igrf([decimal_year,0,float(row['inter_lat']),float(row['inter_lon'])])
        alpha = float(row['strike']) - igrf[0]
        e = np.degrees(np.arctan2(np.tan(np.deg2rad(igrf[1])),np.sin(np.deg2rad(alpha))))
        aei = wrap_180_180(180 - e - wrap_180_180(float(row['phase_shift'])) + asf(srf(row['sz_name'],(float(row['age_max'])+float(row['age_min']))/2)))
        deskew_df.set_value(i,'ei',e)
        deskew_df.set_value(i,'aei',aei)

    return deskew_df

def get_asf_srf(spreading_rate_model_path=None,anomalous_skewness_model_path=None):
    if spreading_rate_model_path==None:
        if os.path.isfile(os.path.join('raw_data','spreading_rate_model.txt')):
            spreading_rate_model_path = os.path.join('raw_data','spreading_rate_model.txt')
        else:
            spreading_rate_model_path = input("spreading rate model file could not be found in raw_data/spreading_rate_model.txt please provide a path to a spreading rate model file here: ")
    if anomalous_skewness_model_path==None:
        if os.path.isfile(os.path.join('raw_data','anomalous_skewness_model.txt')):
            anomalous_skewness_model_path = os.path.join('raw_data','anomalous_skewness_model.txt')
        else:
            anomalous_skewness_model_path = input("anomalous skewness model file could not be found in raw_data/anomalous_skewness_model.txt please provide a path to a anomalous skewness model file here: ")

    srf,sz_list = generate_spreading_rate_model(spreading_rate_model_path)
    asf = generate_anomalous_skewness_model(anomalous_skewness_model_path)

    return asf,srf,sz_list

def phase_shift_data(mag_data,phase_shift):
    mag = np.array(list(map(float,mag_data)))

    N1 = len(mag)

    #   Add a head and a tail for the data
    #   so the magnetic anamaly profile will have two smooth ends
    magpre=np.arange(0,1,0.01)*mag[0]
    magpost=np.arange(.99,-.01,-.01)*mag[-1]
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
    return NewMag[100:N1+100]

def correct_cande(cande_cor_path,deskew_path,spreading_rate_path=os.path.join('raw_data','spreading_rate_model.txt'),dist_e=.75):
    #backup the .deskew file
    if not os.path.isfile(deskew_path+'.ccbak'):
        print("backing up %s to %s"%(deskew_path,deskew_path+'.ccbak'))
        shutil.copyfile(deskew_path,deskew_path+'.ccbak')

    #read in the deskew and cande_cor file
    deskew_df = pd.read_csv(deskew_path,sep="\t")
    cande_cor_df = pd.read_csv(cande_cor_path,sep="\t")
    spreading_rate_func,sz_list = generate_spreading_rate_model(spreading_rate_path)

    #copy the deskew df so I can change it and save the corrected latitudes and longitudes
    new_deskew_df = deskew_df.copy()
    for i,drow in deskew_df.iterrows():
        crow = cande_cor_df[cande_cor_df["comp_name"]==drow["comp_name"]]

        if drow['comp_name'].startswith('#'): continue #commented lines check
        if crow.empty: print("no correction found for component %s"%drow["comp_name"]); continue #no correction for this component check

        half_age = (drow['age_max']-drow['age_min'])/2
        avg_age = (drow['age_max']+drow['age_min'])/2
        half_dis = -np.sign(float(crow['correction']))*half_age*spreading_rate_func(drow['sz_name'],avg_age)

        if drow['track_type'] == 'aero':
            #Find other component direction so we can average the shift between components
            if 'E' in drow["comp_name"]:
                other_crow = cande_cor_df[cande_cor_df["comp_name"]==drow["comp_name"].replace('E','V')]
            elif 'V' in drow["comp_name"]:
                other_crow = cande_cor_df[cande_cor_df["comp_name"]==drow["comp_name"].replace('V','E')]
            else: print("Problem determining component direction for %s"%drow["comp_name"]); continue
            #check that the intersept distance correction between E and V are not more than 3 deg different
            if abs(float(crow['correction']) + float(other_crow['correction']))>3:
                print("correction for %s is >3 km different from the other componenet's correction, and the average may be off"%(drow['comp_name']))
            correction = (float(crow['correction'])+float(other_crow['correction']))/2 + half_dis
        elif drow['track_type'] == 'ship':
            correction = float(crow['correction']) + half_dis
        else:
            print("could not determine the track type for %s please check your deskew file, skipping"%drow["comp_name"]); continue

        corrected_lon,corrected_lat,corrected_dist = get_lon_lat_from_plot_pick(drow,correction,dist_e=dist_e)

        new_deskew_df.set_value(i, 'inter_lat', corrected_lat)
        new_deskew_df.set_value(i, 'inter_lon', corrected_lon)

    new_deskew_df.to_csv(deskew_path,sep="\t",index=False)

def get_lon_lat_from_plot_pick(deskew_row,plot_pick,dist_e=.75):
    drow,correction=deskew_row,plot_pick

    data_file_path = os.path.join(drow["data_dir"],drow["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","idk","mag","lat","lon"],delim_whitespace=True)

    projected_distances = calc_projected_distance(drow['inter_lon'],drow['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),drow['strike'])

    found_dist=False
    for j,row in projected_distances.iterrows():
        if row['dist']>=correction-dist_e and row['dist']<=correction+dist_e:
            picked_lat = round(row['lat'],3) #change lat for the new deskew file
            picked_lon = round(convert_to_0_360(row['lon']),3) #change lon for the new deskew file
            picked_distance = row['dist']
            found_dist=True
            break

    if found_dist:
        print("found lat lon of %s at a distance %.3f"%(drow["comp_name"],picked_distance))
    else:
        print("couldn't find picked distance in datafile to calculate lat and lon for %s"%drow["comp_name"]); return (drow['inter_lon'],drow['inter_lat'],0)

    return picked_lon,picked_lat,picked_distance

def create_maxtab_file(deskew_path,anomoly_name):
    deskew_df = pd.read_csv(deskew_path,sep="\t")
    dates_file = open("raw_data/aeromagCDinfo.txt",'r')
    dates_data = {line.split('\t')[0]:line.split('\t')[3] for line in dates_file.readlines()}
    out_str = ""
    for i,row in deskew_df.iterrows():
        if '#' in row['comp_name']: continue
        track_name = row['comp_name'].split('.')[0]
        if row['track_type']=='ship': date = "%.2f"%get_shipmag_decimal_year(row)
        else: date = "%.2f"%float(dates_data[track_name])
        phase_shift = "%.2f"%convert_to_0_360(float(row['phase_shift']))
        out_str += ' '*(17-len(row['comp_name']))+ row['comp_name'] + ' '
        out_str += 'V' if 'Vd' in row['comp_name'] else 'E'
        out_str += ' '*(6-len(anomoly_name)) + str(anomoly_name)
        out_str += ' '*(9-len(date)) + date
        out_str += ' '*(8-len("%.2f"%float(row['inter_lat']))) + "%.2f"%float(row['inter_lat'])
        out_str += ' '*(8-len("%.2f"%float(row['inter_lon']))) + "%.2f"%float(row['inter_lon'])
        out_str += ' '*(8-len("%.2f"%(float(row['strike'])))) + "%.2f"%(float(row['strike']))
        out_str += ' '*(7-len(phase_shift)) + phase_shift
        out_str += ' '*(11-len('10.000')) + '10.000' + '\n'
    print("saving to maxtab.%s"%anomoly_name)
    out_file = open("maxtab.%s"%anomoly_name,'w+')
    out_file.write(out_str)
    out_file.close()

def generate_synth(ship_or_aero,age_min,age_max,spreading_rate_path,age_path=None,synth_age_lb=-83.64,synth_age_ub=83.64):
    if ship_or_aero=='aero':
        synth_age_path = os.path.join('raw_data','synout_aero')
    elif ship_or_aero=='ship':
        synth_age_path = os.path.join('raw_data','synout_ship')
    elif os.path.isfile(ship_or_aero):
        synth_age_path = ship_or_aero; ship_or_aero = 'ship' #default calls it a ship synthetic
    else: print("generate synth wasn't told which synthetic to generate ship or aero default or path to raw data file, was given %s instead"%str(ship_or_aero)); return

    synth_age = np.loadtxt(synth_age_path)
    NDIM = len(synth_age) #use buffer calculation from synthanom
    ilb = int(((NDIM-1)*(.15/1.3)) + 1) - 1 #minus 1 because python indexes from 0 unlike fortran :p
    iub = int((NDIM-1)*(1.15/1.3) + 0.5) - 1
    index_syn=synth_age[ilb:iub,0]
    mag_syn=synth_age[ilb:iub,1]*1e5 #converting from Gauss to nT (idk why, CGS master race!!!!)

    age_syn=((synth_age_lb-synth_age_ub)/(iub-ilb))*(index_syn-(ilb+1))+synth_age_ub

    spreading_rate_func,sz_list = generate_spreading_rate_model(spreading_rate_path)

    for sz in sz_list:
        print("generating %smag synthetic for spreading zone %s"%(ship_or_aero,sz))
        sz_spreading_rate_func = lambda age: spreading_rate_func(sz,age)
        dis_syn=np.zeros(len(age_syn))
        for num in range(1,len(age_syn)):
            dis_syn[num]=dis_syn[num-1]+(age_syn[num-1]-age_syn[num])*sz_spreading_rate_func(age_syn[num])

        step=-.0001

        sr_max_range = np.array(list(map(sz_spreading_rate_func,np.arange(max(age_syn),age_max+step,step))))
        dis_max = sum(abs(step)*sr_max_range)

        sr_min_range = np.array(list(map(sz_spreading_rate_func,np.arange(age_max,age_min+step,step))))
        dis_min = dis_max+sum(abs(step)*sr_min_range)

        center_dis=(dis_min+dis_max)/2
        dis_min=dis_min-center_dis
        dis_max=dis_max-center_dis
        dis_syn=dis_syn-center_dis
        dis_span=np.array([dis_min,dis_max])

        check_dir('SynthData')
        np.savetxt(os.path.join('SynthData','dis_syn_%s_%s.txt'%(sz,ship_or_aero)),dis_syn,delimiter=',')
        np.savetxt(os.path.join('SynthData','mag_syn_%s_%s.txt'%(sz,ship_or_aero)),mag_syn,delimiter=',')
        np.savetxt(os.path.join('SynthData','dis_span_%s.txt'%sz),dis_span,delimiter=',')

        if age_path==None: age_path=os.path.join('raw_data','timescale_gradstein2012.txt')
        if os.path.isfile(age_path):
            ages_df = pd.read_csv(age_path,sep='\t',header=0,index_col=0)

            step = .0001
            prev_anom = ages_df.index[0]
            ages_df.set_value(ages_df.index[0],'age_span',ages_df['base'][prev_anom])
            ages_df.set_value(ages_df.index[0],'dist',0)
            for anom,age in ages_df['base'].iloc[1:].iteritems():
                age_span = age-ages_df['base'][prev_anom]
                ages_df.set_value(anom,'age_span',age_span)
                age_span_array = np.arange(ages_df['base'][prev_anom],age+step,step)
                sr_array = np.array(list(map(sz_spreading_rate_func,age_span_array)))
                dist = ages_df['dist'][prev_anom]+sum(abs(step)*sr_array)
                ages_df.set_value(anom,'dist',dist)
                prev_anom = anom

#        prev_anom = ages_df.index[0]
#        ages_df.set_value(ages_df.index[0],'age_span',ages_df['base'][prev_anom])
#        for anom,age in ages_df['base'].iloc[1:].iteritems():
#            age_span = age-ages_df['base'][prev_anom]
#            ages_df.set_value(anom,'age_span',age_span)
#            prev_anom = anom

#        prev_anom = ages_df.index[-1]
#        ages_df.set_value(ages_df.index[-1],'dist',0)
#        for anom,age_span in ages_df['age_span'].iloc[-2::-1].iteritems():
#            dist = ages_df['dist'][prev_anom]+ages_df['age_span'][prev_anom]*spreading_rate_func(ages_df['base'][prev_anom])
#            ages_df.set_value(anom,'dist',dist)
#            prev_anom = anom

        [ages_df.set_value(anom,'dist',(max(ages_df['dist'].tolist())-dist)-center_dis) for anom,dist in ages_df['dist'].iteritems()]

        ages_df['dist'].to_csv(os.path.join('SynthData','anom_spans_%s.txt'%sz),sep=',')

def check_generate_synth(sz_to_check,age_min,age_max,ship_or_aero="both",matlab=False,spreading_rate_path=None,synth_data_location='.'):

    if ship_or_aero!="aero" and ship_or_aero!="both" and ship_or_aero!="ship":
        raise ValueError("synthetic type must be either both, aero, or ship, but was given %s please check your inputs"%str(ship_or_aero))

    if ship_or_aero=="aero" or ship_or_aero=="both":

        if not os.path.isfile(os.path.join(synth_data_location,'SynthData/dis_span_%s.txt'%sz_to_check)) or \
           not os.path.isfile(os.path.join(synth_data_location,'SynthData/dis_syn_%s_aero.txt'%sz_to_check)) or \
           not os.path.isfile(os.path.join(synth_data_location,'SynthData/mag_syn_%s_aero.txt'%sz_to_check)):
            if spreading_rate_path==None:
                if os.path.isfile(os.path.join('raw_data','spreading_rate_model.txt')):
                    spreading_rate_path = os.path.join('raw_data','spreading_rate_model.txt')
                else:
                    spreading_rate_path = input("couldn't find spreading rate model file, please provide a path here: ")
            generate_synth('aero',float(age_min),float(age_max),spreading_rate_path)
        else: print("Synthetic files found skipping synthetic generation to save time, please move/remove these files if you would like to force this script to regenerate the synthetic")

    if ship_or_aero=="ship" or ship_or_aero=="both":

        if not os.path.isfile(os.path.join(synth_data_location,'SynthData/dis_span_%s.txt'%sz_to_check)) or \
           not os.path.isfile(os.path.join(synth_data_location,'SynthData/dis_syn_%s_ship.txt'%sz_to_check)) or \
           not os.path.isfile(os.path.join(synth_data_location,'SynthData/mag_syn_%s_ship.txt'%sz_to_check)):
            if spreading_rate_path==None:
                if os.path.isfile(os.path.join('raw_data','spreading_rate_model.txt')):
                    spreading_rate_path = os.path.join('raw_data','spreading_rate_model.txt')
                else:
                    spreading_rate_path = input("couldn't find spreading rate model file, please provide a path here: ")
            generate_synth('ship',float(age_min),float(age_max),spreading_rate_path)
        else: print("Synthetic files found skipping synthetic generation to save time, please move/remove these files if you would like to force this script to regenerate the synthetic")

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
        for ub_sr,as_val in as_lists:
            if float(sr)<=ub_sr: return as_val
        print("spreading rate %.3f could not be found in the anomalous skewness model"%float(sr))

    return anomalous_skewness_model

def get_shipmag_decimal_year(row,deg_e=.01):
    """
    takes a row (pandas series) of a deskew file which is of type ship and returns the decimal year for the
    intersection point. returns none if not found.
    """
    if row['track_type']!='ship':
        raise ValueError("get_shipmag_decimal_year can only run on shipmag data recieved data of type %s instead"%str(row['track_type']))
    data_file_path = os.path.join(row["data_dir"],row["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","decimal_year","mag","lat","lon"],delim_whitespace=True)
    decimal_year=None
    for j,datarow in data_df.iterrows(): #iterate to find the distance associated with the current lat lon
        if (float(datarow['lat'])>=float(row['inter_lat'])-deg_e and \
          float(datarow['lat'])<=float(row['inter_lat'])+deg_e) and \
          (convert_to_0_360(datarow['lon'])>=convert_to_0_360(row['inter_lon'])-deg_e and \
          convert_to_0_360(datarow['lon'])<=convert_to_0_360(row['inter_lon'])+deg_e):
            decimal_year=float(datarow['decimal_year']); break
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
        next_sum_sq = (float(row['inter_lat'])-float(datarow['lat']))**2 + (convert_to_0_360(row['inter_lon'])-convert_to_0_360(datarow['lon']))**2
        if sum_sq<next_sum_sq: break
        elif (float(datarow['lat'])>=float(row['inter_lat'])-deg_e and \
          float(datarow['lat'])<=float(row['inter_lat'])+deg_e) and \
          (convert_to_0_360(datarow['lon'])>=convert_to_0_360(row['inter_lon'])-deg_e and \
          convert_to_0_360(datarow['lon'])<=convert_to_0_360(row['inter_lon'])+deg_e):
            sum_sq = next_sum_sq
            decimal_year=float(datarow['decimal_year'])
    return decimal_year

def reduce_to_pole(deskew_path, pole_lon, pole_lat, spreading_rate_model_path=None, anomalous_skewness_model_path=None):

    asf,srf,sz_list = get_asf_srf(spreading_rate_model_path,anomalous_skewness_model_path)

    deskew_df = filter_deskew_and_calc_aei(deskew_path)

    print("reducing to pole - lat: %.3f, lon: %.3f"%(pole_lat,pole_lon))

    for i,row in deskew_df.iterrows():
        #taken from Lin's Matlab Code (doesn't seem to work)
#        geodes_inv_dic = Geodesic.WGS84.Inverse(float(row['inter_lat']),float(row['inter_lon']),pole_lat,pole_lon)
#        strike,arc_len_pole,az_pole = float(row['strike']),geodes_inv_dic['a12'],geodes_inv_dic['azi1']
#        I = np.arctan(2*np.tan(np.deg2rad(90-arc_len_pole)))
#        a = az_pole+180-strike
#        e = wrap_180_180(np.degrees(np.arctan2(np.tan(np.deg2rad(I)),np.sin(np.deg2rad(a)))))
#        reduced_skewness = wrap_180_180(float(row['phase_shift']) - (e - float(row['ei'])))

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

        reduced_skewness = wrap_0_360(180 - float(row['ei']) - e_r + anom_skew)

        deskew_df.set_value(i,'phase_shift',round(reduced_skewness,3))

    old_results_dir = deskew_df['results_dir'].iloc[0]
    new_results_dir = os.path.join(old_results_dir,"pole_%.0f_%.0f_results"%(pole_lon,pole_lat))
    check_dir(new_results_dir)
    deskew_df['results_dir'] = new_results_dir

    reduced_deskew_df = deskew_df[['comp_name','phase_shift','step','age_min','age_max','inter_lat','inter_lon','strike','data_dir','results_dir','track_type','sz_name','r','g','b']]

    out_path = os.path.join(os.path.dirname(deskew_path),"pole_%.0f_%.0f.deskew"%(pole_lon,pole_lat))   
    print("writing to %s"%out_path)
    reduced_deskew_df.to_csv(out_path,sep='\t',index=False)

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
    colors,i,j = [(1,0,0),(0,1,0),(0,0,1),(.58,0,.83),(0,0,0),(0,1,1),(.5,.5,.5),(1,1,0),(.8,0,1),(.5,0,0),(0,.5,0),(0,0,.5),(.5,0,.5),(0,.5,.5),(.5,.5,0)],0,1
    sz_to_color,sz_to_name = {},{}
    for track,sz,inter in track_sz_and_inters:

        idx = read_idx_from_string(inter)
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
        azsz_df = pd.read_csv(azsz_path,sep=' ')
        strike = azsz_df['strike'][0]

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

    check_dir(data_directory)
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
    deskew_df = pd.read_csv(deskew_path,sep='\t')
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
        check_dir(out_dir)
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
        dists = [Geodesic.WGS84.Inverse(ar[0,1],convert_to_0_360(ar[0,0]),lat,convert_to_0_360(lon))['s12']/1000 for lat,lon in zip(ar[:,1],convert_to_0_360(ar[:,0]))]
        new_dists = np.arange(0,dists[-1],10)
        ridge_lons = np.interp(new_dists,dists,convert_to_0_360(ar[:,0]))
        ridge_min_lon = min(ridge_lons)
        ridge_max_lon = max(ridge_lons)
        ridge_lats = np.interp(new_dists,dists,ar[:,1])
        ridge_min_lat = min(ridge_lats)
        ridge_max_lat = max(ridge_lats)
        ridge_lon_lats = [[lon,lat] for lon,lat in zip(ridge_lons,ridge_lats)]
        ridge_loc_dict[sz] = ridge_lon_lats

    def get_ridge_loc(sz,track_lon_lats):
        if sz not in ridge_loc_dict.keys(): print("sz not found when looking for ridge location, was given spreading zone %s but only had options %s"%(str(sz),str(ridge_loc_dict.keys()))); return None,None
        idx = intersect_bf(track_lon_lats, ridge_loc_dict[sz],e=.5)
        if idx == [None,None]: print("could not calculate intersept"); return None,None
        else: return idx[0][0]

    return get_ridge_loc

def read_and_fit_fz_data(fz_directory=os.path.join('raw_data','fracture_zones')):
    fzs = glob.glob(os.path.join(fz_directory,'*'))
    lfz = []
    for fz in fzs:
        fzdf = pd.read_csv(fz,sep='\t')
        dists = [Geodesic.WGS84.Inverse(fzdf['Latitude'][0],convert_to_0_360(fzdf['Longitude'][0]),lat,convert_to_0_360(lon))['s12']/1000 for lat,lon in zip(fzdf['Latitude'],fzdf['Longitude'])]
        new_dists = np.arange(0,dists[-1],10)
        fz_lons = np.interp(new_dists,dists,convert_to_0_360(fzdf['Longitude']))
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

def find_fz_crossings(deskew_path,fz_directory=os.path.join('raw_data','fracture_zones')):
    deskew_df = pd.read_csv(deskew_path,sep='\t')
    get_fz_loc = read_and_fit_fz_data(fz_directory)

    fz_inter_dict = {}
    for i,row in deskew_df.iterrows():
        data_path = os.path.join(row['data_dir'],row['comp_name'])
        data_df = open_mag_file(data_path)
        track_lon_lats = [[convert_to_0_360(lon),lat] for lon,lat in zip(data_df['lon'],data_df['lat'])]
        inters = get_fz_loc(track_lon_lats)
        if inters != []: fz_inter_dict[row['comp_name']] = inters

    fz_inter_df = pd.DataFrame({'inters':fz_inter_dict})
    fz_inter_df.to_csv('fz_intercepts.txt',sep='\t')

def update_useable_tracks_from_deskew(deskew_path,useable_track_path):
    useable_df = pd.read_csv(useable_track_path, sep='\t', header=None)
    useable_df['tracks'] = list(map(os.path.basename, useable_df[0].tolist()))
    deskew_df = pd.read_csv(deskew_path,sep='\t')
    useable_tracks = list(map(lambda x: x.rstrip('.Ed .lp .Vd'), deskew_df['comp_name'].tolist()))
    new_useable_df = useable_df[useable_df['tracks'].isin(useable_tracks)][[0,1,2]]
    directory = os.path.dirname(useable_track_path)
    new_useable_track_filename = 'new_' + os.path.basename(useable_track_path)
    out_path = os.path.join(directory,new_useable_track_filename)
    new_useable_df.to_csv(out_path, sep='\t', index=False, header=False)

def create_deskewed_data_file(deskew_path):
    #read deskew file
    deskew_df = filter_deskew_and_calc_aei(deskew_path)

    #iterate mag files
    for i,row in deskew_df.iterrows():
        #read mag files
        data_path = os.path.join(row['data_dir'],row['comp_name'])
        data_df = open_mag_file(data_path)
        #deskew mag data
        data_df['deskewed_mag'] = phase_shift_data(data_df['mag'],float(row['phase_shift']))
        #save deskewed mag data as $DATAFILE.deskew
        print("writing %s"%(data_path+'.deskewed'))
        data_df[['lon','lat','deskewed_mag']].to_csv(data_path+'.deskewed',sep=',',header=False,index=False)

