import os
import pandas as pd
import numpy as np
from geographiclib.geodesic import Geodesic

def check_dir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def remove_off_map_points(l):
    return list(filter(lambda x: x<1e30 and x>-1e30, l))

def convert_to_0_360(l):
    cor = lambda x: 360+float(x) if float(x)<0 else float(x)
    if isinstance(l,str) or isinstance(l,int) or isinstance(l,float): return cor(l)
    else: return list(map(cor, l))

def convert_to_180_180(l):
    cor = lambda x: float(x)-360 if float(x)>180 else float(x)
    if isinstance(l,str) or isinstance(l,int) or isinstance(l,float): return cor(l)
    else: return list(map(cor, l))

def wrap_180_180(x):
    try: x = float(x)
    except ValueError: raise ValueError("could not coerce input %s into type float for calculation"%str(x))
    if x>180: return -(-x%180)
    elif x<-180: return (x%180)
    else: return x

def wrap_0_360(x):
    try: x = float(x)
    except ValueError: raise ValueError("could not coerce input %s into type float for calculation"%str(x))
    return x%360

def read_idx_from_string(idx_str):
    flat_idx = list(map(float,idx_str.replace("[","").replace("(","").replace("]","").replace(")","").replace("'","").split(",")))
    return (((flat_idx[0],flat_idx[1]),(flat_idx[2],flat_idx[3])),(int(flat_idx[4]),int(flat_idx[5])))

def calc_projected_distance(inter_lon,inter_lat,lons,lats,strike):
    geodesic_solutions = [Geodesic.WGS84.Inverse(float(inter_lat),convert_to_0_360(inter_lon),lat,convert_to_0_360(lon)) for lat,lon in zip(lats,lons)]

    projected_distances = pd.DataFrame([{'lat':gs['lat2'],'lon':gs['lon2'],'dist':(gs['s12']*np.sin(np.deg2rad(float(strike)-gs['azi2'])))/1000} for gs in geodesic_solutions])

    return projected_distances

def open_deskew_file(deskew_path):
    deskew_df = pd.read_csv(deskew_path,sep='\t')
    deskew_df = deskew_df[~deskew_df["comp_name"].str.startswith('#')]
    cols = deskew_df.columns
    return deskew_df.reset_index()[cols]

def open_mag_file(mag_file):
    dfin = open_aeromag_file(mag_file)
    if dfin.empty: dfin = open_shipmag_file(mag_file)
    if dfin.empty: print("could not open %s as either aeromag or shipmag file, check your data"%mag_file)
    return dfin

def open_shipmag_file(shipmag_file):
    fin = open(shipmag_file,'r')
    lines = fin.readlines()
    fin.close()
    lines = [line.split() for line in lines]
    try: dfin = pd.DataFrame(lines,columns=["dist","dec_year","mag","lat","lon"])
    except AssertionError:
#        print("ship mag file %s does not have the standard 5 rows, you should check this data and see if something happened during processing. Returning a empty dataframe"%shipmag_file)
        dfin = pd.DataFrame()
    return dfin

def open_aeromag_file(aeromag_file):
    fin = open(aeromag_file,'r')
    lines = fin.readlines()
    fin.close()
    lines = [line.split() for line in lines]
    try: dfin = pd.DataFrame(lines,columns=["time","lat","lon","n_comp","e_comp","h_comp","v_comp","mag","dec","inc","None","alt"])
    except AssertionError:
#        print("aeromag file %s does not have the standard 12 rows, you should check this data and see if something happened during processing. Returning a empty dataframe"%aeromag_file)
        dfin = pd.DataFrame()
    return dfin

def write_mag_file_df(df,path):
    df.to_csv(path, sep=' ', header=False, index=False)

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def get_barckhausen_2013_chrons():
    ndf = open('raw_data/chrons/Barckhausen2013/GSFML.Barckhausen++_2013_MGR.picks.gmt','r')

    meta_lonlat = ndf.readlines()[7:]

    m = (9-3)/(-110-(-145))
    b = 3-m*(-145)+1
    ccz_div_line = lambda x: m*x+b #yes I picked this line by hand as it was all I could think of quit judging me. -_-

    ccz_chrons = {}
    gcz_chrons = {}
    for line in meta_lonlat:
        if line.startswith('#'):
            chron = line.lstrip('# @DC').split('|')[0].replace('n','').replace('.1','')
            if chron not in list(ccz_chrons.keys()): ccz_chrons[chron] = []
            if chron not in list(gcz_chrons.keys()): gcz_chrons[chron] = []
            else: continue
        else:
            nd = list(map(float,line.split('\t')))
            if nd[1]>ccz_div_line(nd[0]): ccz_chrons[chron].append(nd)
            else: gcz_chrons[chron].append(nd)

    return ccz_chrons, gcz_chrons

def format_coord(x, y):
    """
    If for some reason the matplotlib backend doesn't visualize your mouse position replace ax.format_coord with 
    this function. According to the internet this should fix the problem anyway.

    Code modifed from user tacaswell on Stack Overflow post:
    https://stackoverflow.com/questions/15876011/add-information-to-matplotlib-navigation-toolbar-status-bar
    """
    col = int(x+0.5)
    row = int(y+0.5)
    return 'x=%1.4f, y=%1.4f'%(x, y)

