import os,sys
import numpy as np
import pandas as pd
import pmagpy.ipmag as ipmag
from pyskew.utilities import check_dir,dt_to_dec,convert_to_0_360
from datetime import datetime
from geographiclib.geodesic import Geodesic


def split_m77t(h77tf,m77tf,data_directory="shipmag_data"):

    check_dir(data_directory)
    h77t_df = pd.read_csv(h77tf,sep='\t',dtype=str)
    m77t_df = pd.read_csv(m77tf,sep='\t',dtype=str)

    new_h77t_files,new_m77t_files = [],[]
    for survey_id in h77t_df['SURVEY_ID']:
        h77t_survey_df = h77t_df[h77t_df['SURVEY_ID']==survey_id]
        m77t_survey_df = m77t_df[m77t_df['SURVEY_ID']==survey_id]

        if "MAG_TOT" not in m77t_survey_df.columns: print("magnetic data not found for %s, skipping"%survey_id); continue
        m77t_survey_df = m77t_survey_df[m77t_survey_df["MAG_TOT"].notnull()] #get only data with uncorrected mag
        m77t_survey_df = m77t_survey_df[m77t_survey_df["DATE"].notnull()] #get only data with dates so they can be corrected

        if m77t_survey_df.empty: print("no magnetic data found in %s, skipping"%survey_id); continue

        survey_dir = os.path.join(data_directory,survey_id)
        check_dir(survey_dir)
        h77t_survey_df.to_csv(os.path.join(survey_dir,survey_id+'.h77t'),sep='\t',index=False)
        m77t_survey_df.to_csv(os.path.join(survey_dir,survey_id+'.m77t'),sep='\t',index=False)
        new_h77t_files.append(os.path.join(survey_dir,survey_id+'.h77t'))
        new_m77t_files.append(os.path.join(survey_dir,survey_id+'.m77t'))

    return new_h77t_files,new_m77t_files

def preprocess_m77t(m77tf,data_directory="shipmag_data"):

    #read in data and initialize empty columns and place holder variables
    m77t_df = pd.read_csv(m77tf,sep='\t',dtype=str)
    m77t_df['DECIMAL_YEAR'] = np.nan
    m77t_df['DIS'] = np.nan
    m77t_df['MAG_COR'] = np.nan
    current_dis,prev_lat_lon = 0,[]

    #check for .lp files existance and skip if it has already been made
    fout_name = os.path.join(data_directory, m77t_df['SURVEY_ID'].iloc[0], m77t_df['SURVEY_ID'].iloc[0]+'.lp')
    if os.path.isfile(fout_name): print(".lp file found for %s, skipping to save time. If you would like to regenerate these files please remove them then rerun the script"%str(m77t_df['SURVEY_ID'].iloc[0])); return

    for i,row in m77t_df.iterrows():

        #create decimal year from datetime
        date = str(row['DATE'])
        if date==str(np.nan): print("no date info for record %d of survey %s, skipping this record"%(i,row['SURVEY_ID'])); continue
        dt_row = datetime(int(date[0:4]),int(date[4:6]),int(date[6:8]))
        dec_year = dt_to_dec(dt_row)
        m77t_df.set_value(i,'DECIMAL_YEAR',round(dec_year,5))

        #calculate distance from last point and add to total distance
        if prev_lat_lon!=[]:
            #/1000 to convert m to km
            current_dis += Geodesic.WGS84.Inverse(float(row['LAT']),float(row['LON']),prev_lat_lon[0],prev_lat_lon[1])['s12']/1000
        prev_lat_lon = [float(row['LAT']),float(row['LON'])]
        m77t_df.set_value(i,'DIS',round(current_dis,5))

        #determine IGRF and remove from uncorrected intensity
        igrf_cor = ipmag.igrf([dec_year,0,float(row['LAT']),float(row['LON'])])[2]
        mag_cor = float(row['MAG_TOT']) - igrf_cor
        m77t_df.set_value(i,'MAG_COR',round(mag_cor,5))

    round3_func = lambda x: round(x,3)
    dis_array = list(map(round3_func,np.arange(float(m77t_df['DIS'].tolist()[0]),float(m77t_df['DIS'].tolist()[-1]),1))) #spacing of 1 km, because I can
    decimal_year_array = list(map(round3_func,np.interp(dis_array,list(map(float,m77t_df['DIS'])),list(map(float,m77t_df['DECIMAL_YEAR'].tolist())))))
    mag_cor_array = list(map(round3_func,np.interp(dis_array,list(map(float,m77t_df['DIS'])),list(map(float,m77t_df['MAG_COR'])))))
    lat_array = list(map(round3_func,np.interp(dis_array,list(map(float,m77t_df['DIS'])),list(map(float,m77t_df['LAT'])))))
    lon_array = list(map(round3_func,np.interp(dis_array,list(map(float,m77t_df['DIS'])),convert_to_0_360(m77t_df['LON']))))

    interp_df = pd.DataFrame({'dis':dis_array,'decimal_year':decimal_year_array,'mag_cor':mag_cor_array,'lat':lat_array,'lon':lon_array})

#    #check distance
#    interp_df['dis_check'] = np.nan
#    current_dis,prev_lat_lon = 0,[]
#    for i,row in interp_df.iterrows():
#        #calculate distance from last point and add to total distance
#        if prev_lat_lon!=[]:
#            #/1000 to convert m to km
#            current_dis += Geodesic.WGS84.Inverse(float(row['lat']),float(row['lon']),prev_lat_lon[0],prev_lat_lon[1])['s12']/1000
#        prev_lat_lon = [float(row['lat']),float(row['lon'])]
#        interp_df.set_value(i,'dis_check',round(current_dis,5))

    #write to .lp file
    print("saving %s"%fout_name)
    interp_df[['dis','decimal_year','mag_cor','lat','lon']].to_csv(fout_name,sep='\t',index=False,header=False)



if __name__=="__main__":
    data_directory='../raw_data/new_shipmag'
    if len(sys.argv)>=4: data_directory=sys.argv[3]
    print("extracting data from main .h77t and .m77t files")
    _,new_m77t_files = split_m77t(sys.argv[1],sys.argv[2],data_directory=data_directory)
    for new_m77t_file in new_m77t_files:
        print("running on %s"%new_m77t_file)
        preprocess_m77t(new_m77t_file,data_directory=data_directory)

