import os,sys
import pandas as pd
import numpy as np
from scipy import odr
from scipy.linalg import invpascal
import scipy.stats as stats
from datetime import datetime
from geographiclib.geodesic import Geodesic
from multiprocessing import Process
from functools import reduce,cmp_to_key

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

def wrap_90_90(x):
    try: x = float(x)
    except ValueError: raise ValueError("could not coerce input %s into type float for calculation"%str(x))
    if x>90: return 90-x%90
    elif x<-90: return -90-(x%-90)
    else: return x

def wrap_0_360(x):
    try: x = float(x)
    except ValueError: raise ValueError("could not coerce input %s into type float for calculation"%str(x))
    return x%360

def dt_to_dec(dt):
    """Convert a datetime to decimal year."""
    year_start = datetime(dt.year, 1, 1)
    year_end = year_start.replace(year=dt.year+1)
    return dt.year + ((dt - year_start).total_seconds() /  # seconds so far
        float((year_end - year_start).total_seconds()))  # seconds in year

def read_idx_from_string(idx_str):
    flat_idx = list(map(float,idx_str.replace("[","").replace("(","").replace("]","").replace(")","").replace("'","").split(",")))
    return (((flat_idx[0],flat_idx[1]),(flat_idx[2],flat_idx[3])),(int(flat_idx[4]),int(flat_idx[5])))

def calc_projected_distance(inter_lon,inter_lat,lons,lats,strike):
    geodesic_solutions = [Geodesic.WGS84.Inverse(lat,convert_to_0_360(lon),float(inter_lat),convert_to_0_360(inter_lon)) for lat,lon in zip(lats,lons)]

    projected_distances = pd.DataFrame([{'lat':gs["lat1"],'lon':gs["lon1"],'dist':(gs["s12"]*np.sin(-np.deg2rad(float(strike)-gs["azi1"])))/1000} for gs in geodesic_solutions])

    return projected_distances

def open_sr_model_file(sr_path):
    sr_file = open(sr_path,'r')
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
    return sr_dict

def write_sr_model_file(sr_dict,sr_path):
    with open(sr_path,'w+') as fout:
        for key in sr_dict.keys():
            fout.write(key+"\n")
            sr_data = sr_dict[key]
            for sr_datum in sorted(sr_data,key=cmp_to_key(lambda x,y: x[0]-y[0])):
                fout.write("%.3f\t%.3f\n"%(float(sr_datum[0]),float(sr_datum[1])))

def open_ellipse_file(ellipse_path):
    elipse_file = open(ellipse_path,"r")
    lon,lat,az,a,b = list(map(float,elipse_file.read().split()))
    return lon,lat,az,a,b

def open_deskew_file(deskew_path):
    deskew_df = pd.read_csv(deskew_path,sep='\t',index_col=False)
    deskew_df = deskew_df[deskew_df["comp_name"].notnull()] #remove data that have no data file record
    deskew_df['quality'] = np.where(deskew_df["comp_name"].str.startswith('#'), "b", "g") #make quality column
    deskew_df["comp_name"] = deskew_df['comp_name'].apply(lambda x: x.lstrip("#")) #remove hashes on commented data
#    deskew_df = deskew_df[~deskew_df["comp_name"].str.startswith('#')]
#    cols = deskew_df.columns
#    return deskew_df.reset_index()[cols]
    return deskew_df

def write_deskew_file(deskew_path,deskew_df): #need to implement with , float_format="%.3f" for prettyness
    deskew_df["comp_name"] = np.where(deskew_df["quality"]=="b", "#" + deskew_df["comp_name"], deskew_df["comp_name"])
    deskew_df.drop("quality",inplace=True,axis=1)
    deskew_df.to_csv(deskew_path,sep="\t",index=False,float_format="%.3f")

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
    try: dfin = pd.DataFrame(lines,columns=["dist","dec_year","mag","lat","lon"],dtype=float)
    except (ValueError,AssertionError) as e:
#        print("ship mag file %s does not have the standard 5 rows, you should check this data and see if something happened during processing. Returning a empty dataframe"%shipmag_file)
        dfin = pd.DataFrame()
    return dfin

def open_aeromag_file(aeromag_file):
    fin = open(aeromag_file,'r')
    lines = fin.readlines()
    fin.close()
    lines = [line.split() for line in lines]
    try: dfin = pd.DataFrame(lines,columns=["time","lat","lon","n_comp","e_comp","h_comp","v_comp","mag","dec","inc","None","alt"],dtype=float)
    except (ValueError,AssertionError) as e:
#        print("aeromag file %s does not have the standard 12 rows, you should check this data and see if something happened during processing. Returning a empty dataframe"%aeromag_file)
        dfin = pd.DataFrame()
    return dfin

def write_mag_file_df(df,path):
    df.to_csv(path, sep=' ', header=False, index=False, float_format="%.3f")

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

def get_barckhausen_2013_chrons(barckhausen_path=os.path.join('..','raw_data','chrons','Barckhausen2013','GSFML.Barckhausen++_2013_MGR.picks.gmt')):
    ndf = open(barckhausen_path,'r')

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

def run_in_parallel(target,args=[],kwargs={}):
        if 'darwin' in sys.platform:
            return target(*args,**kwargs)
        else:
            process = Process(target = target,args=args,kwargs=kwargs)
            return process.start()

def polyfit(x,y,degree,err=None,full=False):
    """
    Does a polynomial fit of any degree using the numpy polyfit routine, which uses SVD to minimize sum sq error.
    Then calculates estimated standard deviation of all model parameters as per Bevington (1969) ch 8.
    As well as the goodness of fit value R^2 to provide more information on the fit than the default numpy function

    Parameters:
        x - the independent variable to compute the regression with respect to, note that x is assumed to have no 
            varriance according to these equations. If this is not the case then pass in err = sqrt(xerr**2 + 
            yerr**2) as per Bevington ch 6.
        y - the dependent variable calculated from x, which is assumed to contain all error
            degree - the degree of the polynomial to be fit
    Optional Parameters:
        err - the error in the dependent variable, should be a vector of equal length to x and y, if None error is assumed unity
        full - boolean which determines the amount of statistics returned

    Returns:
        if full:
            pols,sds,chi2,deg_free,r2,res,rank,sv,rcond
        else:
            pols,sds,chi2,deg_free,r2
    """
    if isinstance(err,type(None)): err = np.ones(len(x)) #equal weight
    x,y,err,N = np.array(x),np.array(y),np.array(err),len(x)

    #calculate fit and chi2
    pols = np.polyfit(x,y,degree,w=(1/err),full=full)
    if full: pols,res,rank,sv,rcond = pols
    yp = np.polyval(pols,x)
    chi2 = sum(((y-yp)/err)**2)
    deg_free = (N-degree-1)
    reduced_chi2 = chi2/deg_free

    #calculate precision parameter (R^2)
#    mean = sum(y)/len(y)
#    ssres = sum((y-yp)**2)
#    sstot = sum((y-mean)**2)
#    r2 = 1 - ssres/sstot
    r2 = (stats.pearsonr(x,y)[0])**2


    #calculate std of each model parameter

    #construct matrix alpha described in Bev 8-23
    alpha = []
    for j,a1 in enumerate(pols):
        alpha.append([])
        for k,a2 in enumerate(pols):
            alpha[j].append(sum((err**-2)*(x**j)*(x**k)))

    #invert alpha to get error matrix which has the varriance of the ith model parameter on it's diagonal as per Bev 8-28
    error_matrix = np.linalg.inv(alpha)
    if (err==np.ones(len(x))).all():
        ss = sum((y-yp)**2)/deg_free
        sds = np.sqrt(ss*np.diag(error_matrix))[::-1]
    else: sds = np.sqrt(np.diag(error_matrix))[::-1]

    if full: return pols,sds,chi2,deg_free,r2,res,rank,sv,rcond
    else: return pols,sds,chi2,deg_free,r2

def polyerr(pols,sds,x,xerr=None):
    if isinstance(xerr,type(None)): xerr = np.zeros(len(x))
    pols,sds,x,xerr = np.array(pols)[::-1],np.array(sds)[::-1],np.array([e if e!=0 else .01 for e in x]),np.array(xerr)
    return np.sqrt(len(x)*sum([((pols[i]*(x**i))**2) * ( ((i*(x**(i-1))*xerr)**2)/(x**2) + (sds[i]**2)/(pols[i]**2) ) for i in range(len(pols))]))

def polyenv(pols,x,yerr,xerr=None,center=None):
    if center==None: center = sum(x)/len(x) #by default center the error env on the mean of the independent var.
    x = np.array(x)

    new_pols,new_sds = polyrecenter(pols,x,yerr,center)

    return polyerr(new_pols,new_sds,x-center,xerr=xerr)

def polyrecenter(pols,x,yerr,center):
    n,pols,x,yerr = len(pols),np.array(pols)[::-1],np.array(x),np.array(yerr)

    P = invpascal(n,kind='lower').T
    A = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if (j-i) < 0: continue
            A[i,j] = P[i,j]*center**(j-i)
    new_pols = (np.linalg.inv(A) @ pols)[::-1]

    alpha = []
    for j,a1 in enumerate(new_pols):
        alpha.append([])
        for k,a2 in enumerate(new_pols):
            alpha[j].append(sum((yerr**-2)*((x-center)**j)*((x-center)**k)))

    error_matrix = np.linalg.inv(alpha)
    new_sds = list(np.sqrt(np.diag(error_matrix)))[::-1]

    return new_pols, new_sds

def odrfit(x,y,sx,sy,func,start_vector):
    model = odr.Model(func)
    data = odr.RealData(x,y,sx=sx,sy=sy)
    odr = odr.ODR(data,model,beta0=start_vector)
    out = odr.run()
    pols = out.beta
    sds = out.sd_beta

    yp = func(pols,x)
    mean = sum(y)/len(y)
    ssres = sum((y-yp)**2)
    sstot = sum((y-mean)**2)
    r2 = 1 - ssres/sstot

    return pols,sds,r2

def format_polyfit(pols,sds,r2,yname='y',xname='x'):
    return yname + ' = ' + reduce(lambda x,y: x+'+'+y,["(%.2f+-%.2f)*%s^%d"%(pol,sd,xname,len(pols)-n-1) for n,(pol,sd) in enumerate(zip(pols,sds))]) + '; R2 = %.2f'%r2



