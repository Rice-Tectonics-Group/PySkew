import os
import shutil
import rdp
import pandas as pd
import numpy as np
from geographiclib.geodesic import Geodesic
from functools import reduce
import scipy.spatial as spatial
import pyskew.plot_geographic as pg
import pyskew.utilities as utl
import pmagpy.ipmag as ipmag
try: from tqdm import tqdm
except ImportError: tqdm = lambda x: x


def intersect(a1,a2,e=1):
    tree = spatial.cKDTree(a1)
    return tree.query(a2, e)

def intersect_bf(a1,a2,e=1):
    cur_r = 1e9
    result = [None,None]
    for i1 in a1:
        for i2 in a2:
            r = np.sqrt((float(i1[0])-float(i2[0]))**2 + (float(i1[1])-float(i2[1]))**2)
            if r <= e and r < cur_r:
                cur_r=r
                result[0] = (i1,i2)
                result[1] = (a1.index(i1),a2.index(i2))
    return result

def angle(dir):
    """
    Credit to unutbu from Stack Overflow
    (https://stackoverflow.com/questions/14631776/calculate-turning-points-pivot-points-in-trajectory-path)

    Returns the angles between vectors.

    Parameters:
    dir is a 2D-array of shape (N,M) representing N vectors in M-dimensional space.

    The return value is a 1D-array of values of shape (N-1,), with each value
    between 0 and pi.

    0 implies the vectors point in the same direction
    pi/2 implies the vectors are orthogonal
    pi implies the vectors point in opposite directions
    """
    dir2 = dir[1:]
    dir1 = dir[:-1]
    return np.arccos((dir1*dir2).sum(axis=1)/(np.sqrt((dir1**2).sum(axis=1)*(dir2**2).sum(axis=1))))

def find_corners(lat_lon_file,tolerance=1,min_angle=np.pi*0.22):
    """
    Credit to unutbu from Stack Overflow
    (https://stackoverflow.com/questions/14631776/calculate-turning-points-pivot-points-in-trajectory-path)

    Runs the Ramer-Douglas-Peucker algorithm to simplify the lat, lon path stored in the specified lat_lon_file.
    It prints out the number of points in the origional path and the number of points in the simplified path.
    Then calculates the points in the simplified path that constitute a turn greater than min_angle.
    
    Parameters
    ----------
    lat_lon_file - file containing lattitude and longitude in that order in columns seperated by whitespace
    tolerance - the tolerance of the RDP algorithm more tolerance=less points kept (default=1)
    min_angle - angle in radians that should constitute a "significant" change in heading
    
    Returns
    -------
    Tuple - (numpy array number_of_points_in_origionalx2 containing the origional data, 
    number_of_points_in_simplifiedx2 containing the simplified data, and an array containing the locations in
    the simplified data where turns exceeding min_angle are found)
    """
    path = os.path.expanduser(lat_lon_file)
    points = np.genfromtxt(path)
    print("original number of points = %d"%len(points))
    x, y = points.T

    # Use the Ramer-Douglas-Peucker algorithm to simplify the path
    # http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
    # Python implementation: https://github.com/sebleier/RDP/
    simplified = np.array(rdp.rdp(points.tolist(), tolerance))

    print("simplified number of points = %d"%len(simplified))
    sx, sy = simplified.T

    # compute the direction vectors on the simplified curve
    directions = np.diff(simplified, axis=0)
    theta = angle(directions)
    # Select the index of the points with the greatest theta
    # Large theta is associated with greatest change in direction.
    idx = np.where(theta>min_angle)[0]+1
    print("number of significant turns = %d"%len(idx))
    
    return points, simplified, idx

def shipmag_preprocess(shipmag_files):
    for shipmag_file in shipmag_files:
        if os.path.basename(shipmag_file).split('.')[-1].startswith('c'):
            shutil.copyfile(shipmag_file,shipmag_file+'.lp')
        ship_df = utl.open_mag_file(shipmag_file)
        latlon_df = ship_df[['lat','lon']]
        latlon_file = shipmag_file + ".latlon"
        latlon_df.to_csv(latlon_file, sep=' ', index=False, header=False)

def aeromag_preprocess(aeromag_files,date_file=os.path.join('..','raw_data','dates.aeromag')):
    for aeromag_file in aeromag_files: #iterate over all aeromag files

        track,extension = os.path.basename(aeromag_file).split('.') #segment the name into parts
        #read data and make a empty dataframe for output data
        adf = utl.open_mag_file(aeromag_file)
        ddf = pd.read_csv(date_file,sep='\t',index_col=0)
        idf = pd.DataFrame(columns=['dis','v_comp','e_comp','h_comp','t_comp'])

        dis = 0
        decimal_year = float(ddf.loc[track]['decimal_year'])
        prev_lat,prev_lon = None,None
        for i,row in adf.iterrows(): #iterate over rows

            row["lon"] = utl.convert_to_0_360(row["lon"])
            adf.at[i,"lon"] = row["lon"]

            try:
                #check for data gaps
                if (row['lat']==None) or (row['lon']==None) or (row['alt']==None) or (row['mag']==None) or (row['v_comp']==None) or (row['e_comp']==None): continue
                #check for absurd values outside of the domain of the varible (this will capture null values of -99999)
                elif (abs(float(row['lat']))>90) or (abs(utl.convert_to_180_180(row['lon']))>180) or (float(row['alt'])<0) or (abs(float(row['mag']))==99999) or (abs(float(row['v_comp']))==99999) or (abs(float(row['e_comp']))==99999): continue
            except ValueError as e: continue #This implies a value which is not convertable to a float as all of these should be floats this datum must be skipped

            if prev_lat!=None and prev_lon!=None: #calculate distance
                dis += Geodesic.WGS84.Inverse(float(row['lat']),float(row['lon']),prev_lat,prev_lon)['s12']/1000
            adf.set_value(i,'dis',dis)

            #calculate and remove IGRF
            dec,inc,mag = ipmag.igrf([decimal_year,float(row['alt'])*0.3048e-3,float(row['lat']),float(row['lon'])])
            res_v_comp = mag*np.sin(np.deg2rad(inc))
            res_e_comp = mag*np.cos(np.deg2rad(inc))*np.cos(np.deg2rad(dec))
            res_h_comp = mag*np.cos(np.deg2rad(inc))
            res_t_comp = mag

            adf.set_value(i,'res_v_comp',float(row['v_comp'])-res_v_comp)
            adf.set_value(i,'res_e_comp',float(row['e_comp'])-res_e_comp)
            adf.set_value(i,'res_h_comp',float(row['h_comp'])-res_e_comp)
            adf.set_value(i,'res_t_comp',float(row['mag'])-res_t_comp)

            prev_lat,prev_lon = float(row['lat']),float(row['lon'])

        adf = adf[(adf['res_e_comp']<3000) & (adf['res_v_comp']<3000) & (adf['res_h_comp']<3000) & (adf['res_t_comp']<3000)]

        #remove a second order polynomial fromm the magnetic data I don't know why but this is something done
        for col in ['res_e_comp','res_h_comp','res_v_comp','res_t_comp']:
            pols = np.polyfit(adf['dis'].tolist(),adf[col].tolist(),2)
            mag_fit = np.polyval(pols,adf['dis'].tolist())
            adf['cor'+col.lstrip('res')] = np.array(adf[col].tolist()) - mag_fit

        #iterpolate and round data
        round3_func = lambda x: round(x,3)
        adf_dis_list = list(map(float,adf['dis'].tolist()))
        adf_lat_list = list(map(float,adf['lat'].tolist()))
        adf_lon_list = list(map(float,adf['lon'].tolist()))
        adf_v_list = list(map(float,adf['cor_v_comp'].tolist()))
        adf_e_list = list(map(float,adf['cor_e_comp'].tolist()))
        adf_h_list = list(map(float,adf['cor_h_comp'].tolist()))
        adf_t_list = list(map(float,adf['cor_t_comp'].tolist()))
        idf['dis'] = list(map(round3_func,np.arange(float(adf['dis'].tolist()[0]),float(adf['dis'].tolist()[-1]),.1))) #spacing of 1 km, because I can
        idf['lat'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_lat_list)))
        idf['lon'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_lon_list)))
        idf['alt'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,list(map(lambda x: float(x)*.3048, adf['alt'].tolist())))))
        idf['v_comp'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_v_list)))
        idf['e_comp'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_e_list)))
        idf['h_comp'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_h_list)))
        idf['t_comp'] = list(map(round3_func,np.interp(idf['dis'].tolist(),adf_dis_list,adf_t_list)))

#        adf[['dis','alt','cor_v_comp','lat','lon']].to_csv(aeromag_file+'.Vd',index=False,header=False,sep='\t')
#        adf[['dis','alt','cor_e_comp','lat','lon']].to_csv(aeromag_file+'.Ed',index=False,header=False,sep='\t')
#        adf[['dis','alt','cor_t_comp','lat','lon']].to_csv(aeromag_file+'.Td',index=False,header=False,sep='\t')
        idf[['dis','alt','v_comp','lat','lon']].to_csv(aeromag_file+'.Vd.lp',index=False,header=False,sep='\t')
        idf[['dis','alt','e_comp','lat','lon']].to_csv(aeromag_file+'.Ed.lp',index=False,header=False,sep='\t')
        idf[['dis','alt','h_comp','lat','lon']].to_csv(aeromag_file+'.Hd.lp',index=False,header=False,sep='\t')
        idf[['dis','alt','t_comp','lat','lon']].to_csv(aeromag_file+'.Td.lp',index=False,header=False,sep='\t')

        if extension.startswith('c'):
            shutil.copyfile(aeromag_file,aeromag_file+'.lp')

        latlon_df = adf[['lat','lon']]
        latlon_file = aeromag_file + ".latlon"
        latlon_df.to_csv(latlon_file, sep=' ', index=False, header=False)

def find_track_cuts(tracks, chrons_info, results_directory, tolerance=1, min_angle=10,plot=False):

    track_cuts={}
    utl.check_dir(os.path.join(results_directory,"turning_points"))

    for track in tracks:

        #Get Track name
        track_name = os.path.basename(track).split('.')[0]
        print(track_name)

        #Run find_corners to simplify the path and find the "significant" turns in the track
        try: points, simplified, idx = find_corners(track+".latlon",tolerance=tolerance,min_angle=min_angle*(np.pi/180))
        except IOError: print("file not found: %s"%(track+".latlon"));continue

        x, y = points.T #get points of track
        sx, sy = simplified.T #get points of simplified track
        track_cuts[track]=[[syx,sxx] for sxx,syx in zip(sx[idx],sy[idx])] #save lon and lat of "significant" turns

        if plot:
            pg.plot_track_cuts(x,y,sx,sy,idx,chrons_info,track_name,results_directory)

    return track_cuts

def seperate_chron_into_spreading_zones(chron_to_analyse):
    #separate chrons into the different spreading zones
    spreading_zone_files = []
    chron,chron_color = chron_to_analyse
    fchron = open("../raw_data/chrons/cande/cande.%s"%str(chron))
    string = fchron.read()
    fchron.close()
    spreading_zones = string.split('>')
    utl.check_dir('spreading_zones')
    for i,spreading_zone in enumerate(spreading_zones):
        if spreading_zone == '': continue
        headerless_spreading_zone = spreading_zone.split('\n')[1:]
        headerless_spreading_zone_string = reduce(lambda x,y: x+'\n'+y,headerless_spreading_zone)
        fchron_out_path =os.path.join('spreading_zones','chron%s_sz%d.txt'%(chron,i))
        fchron_out = open(fchron_out_path,'w+')
        fchron_out.write(headerless_spreading_zone_string)
        fchron_out.close()
        spreading_zone_files.append(fchron_out_path)
    ccz,gcz = utl.get_barckhausen_2013_chrons()
    if str(chron) in ccz.keys():
        i+=1
        ccz_data = np.array([(utl.convert_to_0_360(lonlat[0]),float(lonlat[1])) for lonlat in ccz[str(chron)]])
        if len(ccz_data)>1:
            nlons = np.arange(min(ccz_data[:,0]),max(ccz_data[:,0]),.025)
            nlats = np.interp(nlons, ccz_data[:,0], ccz_data[:,1])
            out_str = reduce(lambda x,y: x+'\n'+y, [str(nlon)+' '+str(nlat) for nlon,nlat in zip(nlons,nlats)])
            fchron_out_path =os.path.join('spreading_zones','chron%s_sz%d.txt'%(chron,i))
            print("Barckhausen data for CCZ included in %s"%fchron_out_path)
            fchron_out = open(fchron_out_path,'w+')
            fchron_out.write(out_str)
            fchron_out.close()
            spreading_zone_files.append(fchron_out_path)
    if str(chron) in gcz.keys():
        i+=1
        gcz_data = np.array([(utl.convert_to_0_360(lonlat[0]),float(lonlat[1])) for lonlat in gcz[str(chron)]])
        if len(gcz_data)>1:
            nlats = np.arange(min(gcz_data[:,1]),max(gcz_data[:,1]),.05)
            nlons = np.interp(nlats, gcz_data[:,1], gcz_data[:,0])
            out_str = reduce(lambda x,y: x+'\n'+y, [str(nlon)+' '+str(nlat) for nlon,nlat in zip(nlons,nlats)])
            fchron_out_path =os.path.join('spreading_zones','chron%s_sz%d.txt'%(chron,i))
            print("Barckhausen data for GCZ included in %s"%fchron_out_path)
            fchron_out = open(fchron_out_path,'w+')
            fchron_out.write(out_str)
            fchron_out.close()
            spreading_zone_files.append(fchron_out_path)
    return spreading_zone_files

def get_track_intersects(chron_to_analyse, tracks_or_cuts, spreading_zone_files, data_directory='.', bounding_lats=(-90,90), bounding_lons=(0,360),e=1):
    """ This function works in 0-360 longitude because otherwise there would be a discontinuty in the Pacific the region of interest """
    chron,chron_color = chron_to_analyse
    chron_name = "chron%s"%(str(chron))
    bound_check_func = lambda x: bounding_lats[0]<float(x[1]) and bounding_lats[1]>float(x[1]) and bounding_lons[0]<float(x[0]) and bounding_lons[1]>float(x[0])
    intersecting_tracks,out_string = [],""
    for track in tqdm(tracks_or_cuts):
        print(track)
        dft = utl.open_mag_file(track)
        if dft.empty: continue
        lt = [[utl.convert_to_0_360(lon),float(lat)] for lon,lat in zip(dft['lon'],dft['lat'])]
        if not list(filter(bound_check_func,lt)): print("track out of bounds, skipping track"); continue

        for spreading_zone_file in spreading_zone_files:
            lsz = [[line.split()[0],line.split()[1]] for line in open(spreading_zone_file).readlines()]
            if not list(filter(bound_check_func,lsz)): continue
            idx = intersect_bf(lt,lsz,e=e)

            if not any(idx): continue
            else:
                print("-----------intersected in bounds-------------")
                intersecting_tracks.append(track)
                out_string +="%s\t%s\t%s\n"%(track,spreading_zone_file,str(idx))
                break

    print("found %d intersecting tracks"%len(intersecting_tracks))
    utl.check_dir(data_directory)
    fout_name=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
    if os.path.isfile(fout_name): print("backing up %s to %s"%(fout_name,fout_name+'.bak')); shutil.copyfile(fout_name,fout_name+'.bak')
    fout = open(fout_name,'w+')
    print("writing to %s"%fout_name)
    fout.write(out_string)
    fout.close()

    return intersecting_tracks, fout_name

def seperate_data(data_directory, usable_tracks_path):
    usable_tracks_file = open(usable_tracks_path,'r')
    usable_tracks = usable_tracks_file.readlines()
    usable_tracks_file.close()
    tracks_sz_and_inter = [track.split('\t') for track in usable_tracks]

    aeromag_directory = os.path.join(data_directory,'aero')
    shipmag_directory = os.path.join(data_directory,'ship')
    utl.check_dir(data_directory);utl.check_dir(aeromag_directory);utl.check_dir(shipmag_directory)
    new_tracks,out_string = [],""
    for track,sz,inter in tracks_sz_and_inter:
        track_dir,track_file = os.path.split(track)
        if 'hi_alt' in track_dir:
            new_track_dir = os.path.join(aeromag_directory, track_file.split('-')[0], track_file.split('-')[1][:-4])
        elif 'ship' in track_dir:
            new_track_dir = os.path.join(shipmag_directory, os.path.basename(track_file).split('.')[0])
        else:
            print("couldn't identify if the data for track %s was ship or aeromag just sticking it in the data directory, this does rely on a specific directory structure being used which while bad works so if you want to mess with how this is done feel free."%track)
            new_track_dir = data_directory
        utl.check_dir(new_track_dir)
        new_track_file = os.path.join(new_track_dir, track_file)
        shutil.copyfile(track,new_track_file)
        out_string += "%s\t%s\t%s"%(new_track_file,sz,str(inter))
        new_tracks.append(new_track_file)

    #write new file paths to usable_tracks_file
    seperated_tracks_path = usable_tracks_path
    seperated_tracks_file = open(seperated_tracks_path, 'w+')
    seperated_tracks_file.write(out_string)
    seperated_tracks_file.close()

    return new_tracks,seperated_tracks_path

def cut_tracks_and_flip(track_cuts, data_directory, heading="east"):

    cut_tracks,flipped_data=[],[]
    for track,cuts in track_cuts.items():
        print("Starting Track: %s"%track)
        directory,path = os.path.split(track)
        dfin = utl.open_mag_file(track)
#        fin = open(track,'r')
#        lines = fin.readlines()
#        fin.close()
#        lines = [line.split() for line in lines]
#        dfin = pd.DataFrame(lines,columns=["time","lat","lon","n_comp","s_comp","h_comp","v_comp","mag","dec","inc","None","alt"])
        lats=list(map(float,dfin['lat'].tolist()))
        lons=list(map(utl.convert_to_0_360,dfin['lon'].tolist()))
        df_segments=[]
        for cut in cuts:
            try: cut_index = [[lon,lat] for lon,lat in zip(lons,lats)].index(list(cut))
            except ValueError as e: import pdb; pdb.set_trace()
            print("cutting track: %s along index: %d"%(track,cut_index))
            df_segments.append(dfin.loc[:cut_index])
            dfin = dfin.loc[cut_index:]
        df_segments.append(dfin)
        i=1
        for df_segment in df_segments:
            if len(df_segment)==0: continue
            if heading=='east':
                flip_bool = (utl.convert_to_0_360(df_segment['lon'].iloc[0])>utl.convert_to_0_360(df_segment['lon'].iloc[-1])) #is heading easterly
            elif heading=='west':
                flip_bool = (utl.convert_to_0_360(df_segment['lon'].iloc[0])<utl.convert_to_0_360(df_segment['lon'].iloc[-1])) #is heading westerly
            elif heading=='north':
                flip_bool = (df_segment['lat'].iloc[0]>df_segment['lat'].iloc[-1]) #is heading northerly
            elif heading=='south':
                flip_bool = (df_segment['lat'].iloc[0]<df_segment['lat'].iloc[-1]) #is heading southerly
            else: print("the heading provided is not a cardinal direction please rerun with this corrected");return
            if flip_bool:
                print("flipping data for cut: %d track: %s such that path is %serly and thus (hopefully) oldest data first"%(i,path,heading))
                df_segment = df_segment.iloc[::-1]
                flipped_data.append(track.split('.')[0]+'.c%d'%i)
            if not os.path.isdir(os.path.join(directory,'c%d'%i)):
                print("making directory %s"%os.path.join(directory,'c%d'%i))
                utl.check_dir(os.path.join(directory,'c%d'%i))
            segment_path = os.path.join(directory,'c%d'%i,path.split('.')[0] + '.c%d'%i)
            i+=1
            print("writing: %s"%segment_path)
            df_segment.to_csv(segment_path, sep='\t', header=False, index=False)
            cut_tracks.append(segment_path)
    f_flipped = open(os.path.join(data_directory,"flipped_data.txt"),'w+')
    f_flipped.write(reduce(lambda x,y: x+'\n'+y, flipped_data))
    f_flipped.close()
    return cut_tracks,flipped_data

def seperate_E_V(cut_tracks):
    """
    Depriciated: just moves E and V into seperate files for analysis
    """
    for cut_track in cut_tracks:
        cut_track_dir,cut_track_path = os.path.split(cut_track)
        e_dir = os.path.join(cut_track_dir,'E')
        v_dir = os.path.join(cut_track_dir,'V')
        if not os.path.isdir(e_dir):
            utl.check_dir(e_dir)
        if not os.path.isdir(v_dir):
            utl.check_dir(v_dir)
        e_mv_files = (cut_track_path+".E",cut_track_path+".Ed",cut_track_path+".Ed.lp",cut_track_path+".Ed.xyz")
        for e_mv_file in e_mv_files:
            shutil.copyfile(os.path.join(cut_track_dir,e_mv_file),os.path.join(e_dir,e_mv_file))
        v_mv_files = (cut_track_path+".V",cut_track_path+".Vd",cut_track_path+".Vd.lp",cut_track_path+".Vd.xyz")
        for v_mv_file in v_mv_files:
            shutil.copyfile(os.path.join(cut_track_dir,v_mv_file),os.path.join(v_dir,v_mv_file))

def generate_az_strike_files(track_sz_and_inters, chron_to_analyse, heading, results_directory, plot=False):
    """results_directory is just used to find a good place to output plots az and strike data are saved to the data_directory along with normal data files"""

    chron,chron_color = chron_to_analyse
    chron_name = "chron%s"%(str(chron))
    tracks,az_files = [],[]
    for track,spreading_zone_file,inter in track_sz_and_inters:
        idx = read_idx_from_string(inter)

        gcp_lon,gcp_lat = open(spreading_zone_file[:-3]+'gcp').readlines()[5].split()[0:2] #GMT has varriable ouutput to fitcircle so some of these files will have the N eigin value pole here and some the South this may need to be adapted

        az = utl.convert_to_0_360(Geodesic.WGS84.Inverse(float(idx[0][1][1]),float(idx[0][1][0]),float(gcp_lat),float(gcp_lon))['azi1'])

        if heading=='east' and az>=180: az -= 180
        elif heading=='west' and az<=180: az += 180
        elif heading=='north' and (az<=270 and az>=90): az = utl.wrap_0_360(180+az)
        elif heading=='south'and (az>=270 or az<=90): az = utl.wrap_0_360(180+az)
        elif heading not in ['east','west','north','south']: print("heading is not a cardinal direction was given %s so az cannot be corrected proceeding with %.3f as azimuth of younging for track %s correct manually in *.azszs if this is incorrect"%(heading,float(az),str(track)))

        strike = utl.wrap_0_360(az + 90)

        print("track", os.path.basename(track))
        print("spreading zone", os.path.basename(spreading_zone_file))
        print("azimuth", az, "strike", strike)

        fout_name = track[:-3]+'_'+track[-2:]+'_'+os.path.basename(spreading_zone_file)[:-4]+'.azszs'
        fout = open(fout_name, 'w+')
        out_str = """track chron spreadingzone intercept_track intercept_chron azimuth strike\n"""
        out_str += """%s %s %s %s %s %s %s"""%(track, chron_name, spreading_zone_file, str(idx[1][0]), str(idx[1][1]), az, strike)
        fout.write(out_str)
        fout.close()

        az_files.append(fout_name); tracks.append(track)

        # Deletes the Matlab/Octave file; not needed if Python is used to find strike
        #if os.path.isfile('.tmp_az.txt'): os.remove('.tmp_az.txt')

        if plot:
            pg.plot_az_strike(track,spreading_zone_file,idx,az,strike,chron_color,chron_name,results_directory,fout_name)

    return tracks, az_files

def remove_failed_data(new_useful_tracks_file,data_directory,failed_data_dir="failed_data"):
    #reads in new cuts or tracks that should be kept
    nfin = open(new_useful_tracks_file,'r')
    nfin_lines = nfin.readlines()
    nfin.close()
    ncts = [nfin_line.split('\t')[0] for nfin_line in nfin_lines]

    utl.check_dir(failed_data_dir)

    for d in [x[0] for x in os.walk(data_directory)]:
        if d == data_directory: continue
        elif not os.path.isdir(d): continue
        elif '/E' in d or '/V' in d or '/E_backup' in d or '/V_backup' in d or '/example_E' in d or '/example_V' in d : continue
        elif all([d not in nct for nct in ncts]):
            print("moving directory %s because it is no longer being analysed"%d)
            shutil.move(d,os.path.join(failed_data_dir,d))

    for nct in ncts:
        if not os.path.isfile(nct): print("ERROR: there was a problem cleaning out failed data and %s was moved despite it still being used in analysis please manually put it back")

