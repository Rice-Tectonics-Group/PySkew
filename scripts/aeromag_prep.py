import os,sys
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
from pyskew.utilities import check_dir,dt_to_dec
from datetime import datetime
from geographiclib.geodesic import Geodesic

def split_m88t(h88tf,m88tf,data_directory="aeromag_data"):
    check_dir(data_directory)
    h88t_df = pd.read_csv(h88tf,sep='\t',dtype=str)
    m88t_df = pd.read_csv(m88tf,sep='\t',dtype=str)

    data_type = ""
    if "X" in h88t_df["PARAMS_CO"].iloc[0] and "Y" in h88t_df["PARAMS_CO"].iloc[0] and "Z" in h88t_df["PARAMS_CO"].iloc[0]:
        data_type = "vector"
    elif "T" in h88t_df["PARAMS_CO"].iloc[0]:
        data_type = "total"
    else:
        raise TypeError("Could not identify the type of aeromagnetic data being prepared for analysis from the h88t PARAMS_CO header, please specify if XYZ vector data exist or if only T (total) component data exists here and try again")

    new_h88t_files,new_m88t_files,dates_list = [],[],[]
    for survey_id in h88t_df['SURVEY_ID']:
        h_survey_df = h88t_df[h88t_df['SURVEY_ID']==survey_id]
        survey_df = m88t_df[m88t_df['SURVEY_ID']==survey_id]

        survey_name = survey_id.split("WORLD")[-1].strip("_ -")
        survey_dir = os.path.join(data_directory,survey_name)
        check_dir(survey_dir)

        if len(survey_df['LINEID'].drop_duplicates()) == 1: raise RuntimeError("This method has yet to be implemented for surveys which do not have any line ids and need to be split based on changes in flight base which should be possible by calculating the distance between all points and looking for jumps in the delta distance value")

        for e in survey_df['LINEID'].drop_duplicates():
            track_name = survey_name+'_'+e
            line_df = survey_df[survey_df['LINEID']==e]

            #get avg decimal year for this flight
            dates_dict = {}
            dates = line_df[line_df['DATE'].notnull()]['DATE'].tolist()
            avg_date = sum(map(lambda x: dt_to_dec(datetime(int(x[0:4]),int(x[4:6]),int(x[6:8]))),dates))/len(dates)
            year = dates[0][0:4]
            month = dates[0][4:6]
            day = dates[0][6:8]
            dates_dict['profile'] = track_name
            dates_dict['decimal_year'] = avg_date
            dates_dict['year'] = year
            dates_dict['month'] = month
            dates_dict['day'] = day
            dates_dict['mean'] = 'mean'
            dates_list.append(dates_dict)

            if data_type=="vector":
                line_df = line_df[(line_df["MAG_X_NRTH"].notnull()) & (line_df["MAG_Y_EAST"].notnull()) & (line_df["MAG_Z_VERT"].notnull())]
            else:
                line_df = line_df[(line_df["MAG_TOTOBS"].notnull())]

            if np.isnan(e): e = "nan"

            line_f_name = os.path.join(survey_dir,track_name+".m88t")
            line_df.to_csv(line_f_name,sep='\t',index=False)
            new_m88t_file.append(line_f_name)

        h_survey_name = os.path.join(survey_dir,survey_id+".h88t")
        h_survey_df.to_csv(h_survey_name,sep='\t',index=False)
        new_h88t_files.append(h_survey_name)

    dates_df = pd.DataFrame(dates_list)
    dates_df.reindex('profile')
    dates_df.to_csv(os.path.join(survey_dir,survey_name+".dates"),sep='\t')

    return new_h88t_files,new_m88t_files,date_files

def preprocess_m88t(m88tf,data_directory="aeromag_data"):
    m88t_df = pd.read_csv(m88tf,sep='\t',dtype=str)

    alt_id = ""
    if "ALT_BAROM" in m88t_df.columns and m88t_df[m88t_df["ALT_BAROM"].notnull()].empty:
        alt_id = "ALT_BAROM"
    elif "ALT_RADAR" in m88t_df.columns and m88t_df[m88t_df["ALT_RADAR"].notnull()].empty:
        alt_id = "ALT_RADAR"
    elif "ALT_GPS" in m88t_df.columns and m88t_df[m88t_df["ALT_GPS"].notnull()].empty:
        alt_id = "ALT_GPS"
    else: raise ValueError("No altitude found in %s from ALT_BAROM,ALT_GPS,ALT_RADAR so IGRF cannot be removed add ALT data somehow")

    m88t_df = m88t_df[m88t_df[alt_id].notnull()]
    m88t_df["ALT"] = map(lambda x: float(x)*3.28084, m88t_df[alt_id])

    dat_df = m88t_df[["TIME", "LAT", "LON", "MAG_X_NRTH", "MAG_Y_EAST", "MAG_HORIZ", "MAG_Z_VERT", "MAG_TOTOBS", "MAG_DECLIN", "MAG_INCLIN", "MAG_TOTCOR", "ALT"]]

    dat_df.fillna(value=-99999)

    dat_df.to_csv(m88tf.split(".")[0]+".DAT",sep='\t',index=False,header=False)

if __name__=="__main__":
    data_directory='./raw_data/new_aeromag'
    if len(sys.argv)>=4: data_directory=sys.argv[3]
    print("extracting data from main .h88t and .m88t files")
    _,new_m88t_files,_ = split_m88t(sys.argv[1],sys.argv[2],data_directory=data_directory)
    for new_m88t_file in new_m88t_files:
        print("running on %s"%new_m88t_file)
        preprocess_m88t(new_m88t_file,data_directory=data_directory)

