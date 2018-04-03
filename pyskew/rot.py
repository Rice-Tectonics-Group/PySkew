from numpy import *
import pandas as pd
import pyskew.plot_skewness as psk
from geographiclib.geodesic import Geodesic
from functools import reduce

################################################PlateCircuit Obect#################################################

class PlateCircuit(object):
    """
    
    """

    def __init__(self):
        pass

#############################################PlateReconstruction Obect#############################################

class PlateReconstruction(object):
    """
    
    """

    def __init__(self,mov_plate,fix_plate,rots,comment=None,sources={}):
        self.mov_plate = str(mov_plate)
        self.fix_plate = str(fix_plate)
        if comment==None: comment = "Motion of %s relative to %s"%(mov_plate,fix_plate)
        self.comment = str(comment)
        self.sources = sources
        if any(list(map(lambda x: not isinstance(x,Rot), rots))):
            raise ValueError("All values in rotsn argument must be Rot objects, was given %s"%str(rots))
        try: self.rots = [rot.backwards_rot() for rot in rots]
        except ValueError: raise ValueError("rots argument must be list-like, was given %s"%str(rots))
        self.rots.sort() #insure self.rots is always sorted by age
        self.calculate_stage_and_reconst_rots()

    #################I/O#################

    @classmethod
    def read_csv(self,fin,sep='\t'):
        with open(fin,'r') as f:
            mov_plate,fix_plate,comment = f.readline().rstrip('\n').split('\t')[:3]
        df = pd.read_csv(fin,header=1,sep='\t')
        rots_list = df[['lat','lon','rot','start_age','stop_age']].values.tolist()
        rots = [Rot(*rot) for rot in rots_list]
        sources = df['source'].to_dict()
        new_sources = {source:rot for source,rot in zip(sources,rots)}
        return PlateReconstruction(mov_plate,fix_plate,rots,comment=comment,sources=new_sources)

    def to_csv(self,fout,sep='\t'):
        df = self.to_df()
        with open(fout,'w+') as f:
            f.write("%s\t%s\t"%(self.mov_plate,self.fix_plate)+self.comment+"\n")
        with open(fout, 'a') as f:
            df.to_csv(f, sep=sep, index=False, float_format="%.3f")

    def stages_to_csv(self,fout,sep='\t'):
        df = self.stages_to_df()
        with open(fout,'w+') as f:
            f.write("%s\t%s\t"%(self.mov_plate,self.fix_plate)+self.comment+"\n")
        with open(fout, 'a') as f:
            df.to_csv(f, sep=sep, index=False, float_format="%.3f")

    #################Casting#################

    def to_df(self):
        df = pd.DataFrame([(r.backwards_rot()).to_dict(add_data={'source':self.sources[r]}) for r in self.rots])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','source']]
        return df

    def reconst_to_df(self):
        df = pd.DataFrame([(r).to_dict(add_data={'source':'this study'}) for r in self.reconst_rots])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','source']]
        return df

    def stages_to_df(self):
        df = pd.DataFrame([(r).to_dict(add_data={'source':'this study'}) for r in self.stage_rots])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','source']]
        return df

    def list_ages(self):
        """
        Return non-zero ages for which there are rotations
        """
        ages=[]
        for rot in self.reconst_rots:
            if rot.age_lb not in ages: ages.append(rot.age_lb)
            if rot.age_ub not in ages: ages.append(rot.age_ub)
        if 0.0 in ages: ages.remove(0.0)
        return ages

    def get_age_bounds(self):
        ages = [0.0]+self.list_ages()
        ages.sort()
        return ages[0],ages[-1]

    def get_reconstruction_ages(self):
        """Gets ages at which uninterpolated reconstruction poles exist"""
        return self.reconst_to_df()['stop_age'].tolist()

    #################mpfin Compatibility#################

    @classmethod
    def parse_mpfin(self,mpfin):
        fin = open(mpfin,'r')
        comment = fin.readline().rstrip('\n')
        fix_plate = fin.readline().split()[0]
        lines = fin.readlines()
        current_rots,current_plates,reconstructions = [],[fix_plate,""],[]
        for line in lines:
            values = line.split()
            if values[0]=='zz': plate1='zz'
            else: plate1,plate2 = values[0].split('-')
            if plate1 in current_plates and plate2 in current_plates:
                current_rots.append(Rot(float(values[2]),float(values[3]),float(values[4]),0,float(values[1])))
            else:
                if fix_plate not in current_plates:
                    fix_plate=reconstructions[-1].mov_plate
                if current_plates[0]==fix_plate:
                    fplate,mplate=current_plates
                    current_rots = [-r for r in current_rots]
                elif current_plates[1]==fix_plate: fplate,mplate=current_plates[::-1]
                else: raise RuntimeError("Couldn't find fixed plate for plate pair %s,%s"%(plate1,plate2))
                if current_rots!=[]: reconstructions.append(PlateReconstruction(mplate,fplate,current_rots,comment=comment))
                if plate1=='zz': break
                current_rots=[Rot(float(values[2]),float(values[3]),float(values[4]),0,float(values[1]))]
                current_plates = [plate1,plate2]
        return reconstructions

    @classmethod #needs to be moved to Circuit object
    def do_mpfin(self,mpfin):
        """
        reads .fin file and does the same opperations as mpfin
        """
        raise NotImplementedError("This is a feature in progress for now you can manually manipulate the reconstructions")

    #################Calculation Functions#################

    def calculate_stage_and_reconst_rots(self):
        """
        As of now assumes that lower age bound of previous rotation is same as other rotation so they cancel
        """
        self.stage_rots = [self.rots[0].backwards_rot()]
        self.reconst_rots = [self.rots[0].backwards_rot()]
        if self.rots[0] not in self.sources: self.sources[self.rots[0]] = '' #handle first rotation's sources

        for rot in self.rots[1:]: #construct reconstruction rotations
            if rot not in self.sources: self.sources[rot] = '' #handle sources
            if rot.age_ub in self.get_reconstruction_ages():
                print("Age %.2f already has an uninterpolated pole constained in reconst keeping only new rotation:\n%s\n%s\nDiscarding first rotation"%(rot.age_ub,str(self[rot.age_ub]),str(rot)))
                del self[rot.age_ub]
            prev_age = self.reconst_rots[-1].age_ub
            self.stage_rots.append(self[prev_age:rot.age_lb]+rot)
            self.reconst_rots.append(self[rot.age_lb]+rot)

        return self.stage_rots,self.reconst_rots

    def old_calculate_stage_and_reconst_rots(self):
        """
        As of now assumes that lower age bound of previous rotation is same as other rotation so they cancel
        """
        self.stage_rots = [self.rots[0].backwards_rot()]
        self.reconst_rots = [self.rots[0].backwards_rot()]
        if self.rots[0] not in self.sources: self.sources[self.rots[0]] = '' #handle first rotation's sources
        for rot in self.rots[1:]:
            prev_rot=self.reconst_rots[-1]
            if rot not in self.sources: self.sources[rot] = '' #handle sources
            if rot==prev_rot:
                print("Equivlent rots found:\n%s\n%s\nDiscarding second rotation"%(str(prev_rot),str(rot)))
            else:
                rrot = rot.backwards_rot()
                self.reconst_rots.append(rrot)
                self.stage_rots.append(prev_rot.forward_rot()+rrot)
        return self.stage_rots,self.reconst_rots

    def interpolate(self,other_reconst):
        """
        Interpolate rotations for the set of unique ages between 2 reconstructions

        Parameters:
            other_reconst -> other reconstruction to interpolate with

        Return:
            -> list of rotations which gives total reconstructions for all unique ages between the 2 reconstions
            -> list of rotations which gives total reconstructions for all unique ages between the 2 reconstions
        """
        if not isinstance(other_reconst,PlateReconstruction):
            raise TypeError("can only interpolate between reconstructions, was given %s"%str(other_reconst))
        ages = merge_lists_without_duplicates(self.list_ages(),other_reconst.list_ages())
        return (self[ages],other_reconst[ages])

    #################Magic Functions#################

    ###########Arithmatic###########

    def __add__(self,other_reconst):
        if not isinstance(other_reconst,PlateReconstruction):
            raise TypeError("addition is only defined between reconstructions, was given %s"%str(other_reconst))
        if self.fix_plate==other_reconst.mov_plate:
            my_rots,oth_rots = self.interpolate(other_reconst)
            new_rots = [oth_rot+my_rot for my_rot,oth_rot in zip(my_rots,oth_rots)]
            return PlateReconstruction(self.mov_plate,other_reconst.fix_plate,new_rots)
        elif self.mov_plate==other_reconst.mov_plate:
            return -self+other_reconst
        elif self.fix_plate==other_reconst.fix_plate:
            return self+-other_reconst
        elif self.mov_plate==other_reconst.fix_plate:
            raise RotationError("The fixed plate of the first reconstruction should be in common with the moving plate of the second reconstruction. There are 2 possible reconstructions given your input reverse one of your reconstructions at least so there is a uniqe solution. Was given:\n%s\n%s"%(str(self),str(other_reconst)))
        else:
            raise ValueError("The two reconstructions being added must have at least 1 plate in common. Recived:\nmoving -> fixed\n%s -> %s\n%s -> %s"%(self.mov_plate, self.fix_plate, other_reconst.mov_plate, other_reconst.fix_plate))

    def simple__add__(self,other_reconst):
        if not isinstance(other_reconst,PlateReconstruction):
            raise TypeError("addition is only defined between reconstructions, was given %s"%str(other_reconst))
        if self.fix_plate==other_reconst.mov_plate:
            my_rots,oth_rots = self.interpolate(other_reconst)
            new_rots = [oth_rot+my_rot for my_rot,oth_rot in zip(my_rots,oth_rots)]
            return PlateReconstruction(self.mov_plate,other_reconst.fix_plate,new_rots)
        elif self.mov_plate==other_reconst.mov_plate or self.fix_plate==other_reconst.fix_plate or self.mov_plate==other_reconst.fix_plate:
            raise RotationError("The fixed plate of the first reconstruction must be the same as the moving plate of the second reconstruction. Was given:\n%s\n%s"%(str(self),str(other_reconst)))
        else:
            raise ValueError("The two reconstructions being added must have at least 1 plate in common. Recived:\nmoving -> fixed\n%s -> %s\n%s -> %s"%(self.mov_plate, self.fix_plate, other_reconst.mov_plate, other_reconst.fix_plate))

    def __sub__(self,other_reconst):
        return -self+-other_reconst

    def __neg__(self):
        new_rots = [-r for r in self.rots]
        return PlateReconstruction(self.fix_plate,self.mov_plate,new_rots,comment=self.comment,sources=self.sources)

    def __invert__(self):
        new_rots = [~r for r in self.rots]
        return PlateReconstruction(self.mov_plate,self.fix_plate,new_rots,comment=self.comment,sources=self.sources)

    ###########Indexing & Sequencing###########

    def __len__(self):
        return len(self.rots)

    def __getitem__(self,_slice):
        #check types and parse input
        if isinstance(_slice,slice): start,stop,step = _slice.start,_slice.stop,_slice.step
        elif isinstance(_slice,float) or isinstance(_slice,int): start,stop,step = _slice,None,None
        elif isinstance(_slice,list) or isinstance(_slice,tuple) or isinstance(_slice,ndarray):
            return [self.__getitem__(start) for start in _slice]
        else: raise TypeError("A slice of a rotation is not defined for %s"%str(_slice))

        #check boundedness of problem
        min_age,max_age = self.get_age_bounds()
        if start < min_age:
            raise ValueError("Start age of %.3f is less than youngest reconstruction for %s relative to %s"%(float(start),self.mov_plate,self.fix_plate))
        if start > max_age:
            print("Start age %.3f greater than oldest reconstruction for %s relative to %s, extending final rotation"%(float(start),self.mov_plate,self.fix_plate))
            return self.reconst_rots[-1][start]

        if stop==None and step==None: #Index case
            tot_rot = Rot(0,0,0,0,0)
            for rot in self.stage_rots:
                if start in rot: return tot_rot + rot[start]
                tot_rot += rot
        elif stop!=None and step==None: #Slice without step case
            return self.__getitem__(start).reverse_time() + self.__getitem__(stop)
        else: #Slice with step case
            rots = []
            for value in range(start,stop+step,step):
                if value==0: continue
                tot_rot = Rot(0,0,0,0,0)
                for rot in self.stage_rots:
                    if value in rot: rots.append(tot_rot + rot[value]); break
                    tot_rot += rot
            return rots

        raise RuntimeError("Couldn't find a rotation or a problem with input in __getitem__ for reconstruction %s this is almost certainly a bug, contact a developer at the github page."%(str(self)))

    def __setitem__(self,idx,value):
        if isinstance(value,Rot): self.rots.append(value); self.rots.sort()
        elif isinstance(idx,slice) and len(value)==3:
            self.rots.append(Rot(*value,idx.start,idx.stop))
            self.rots.sort()
        else:
            raise TypeError("Cannot add roation to PlateReconstruction for index: %s and value: %s. You must provide either a Rot object for value and any index (index will be ignored) or a slice for index (i.e. [0:10]) which represents your age bounds and a list with values [lat,lon,w] for your value."%(str(idx),str(value)))
        self.calculate_stage_and_reconst_rots()

    def __delitem__(self,value):
        if isinstance(value,Rot):
            if value in self.rots:
                self.rots.remove(value)
            else: raise ValueError("There is no equivlent rotation in this reconstruction cannot remove:\n%s"%str(value))
        elif isinstance(value,int) or isinstance(value,float):
            for rot in self.rots:
                if value in rot:
                    self.rots.remove(rot)
        else:
            raise TypeError("There is no way to remove %s from reconstruction"%str(value))
        self.calculate_stage_and_reconst_rots()

    def __contains__(self,value):
        if isinstance(value,float) or isinstance(value,int):
            return any([(value in rot) for rot in self.rots])
        elif isinstance(value,Rot):
            return ((value in self.rot) or (value in self.stage_rots) or (value in self.reconst_rots))
        else:
            raise TypeError("cannot check if %s is in a reconstruction object"%str(value))

    ###########Visualizing###########

    def __str__(self):
        min_age,max_age = self.get_age_bounds()
        return "Reconstruction Rotations for %s relative to %s, from %.2f to %.2f"%(self.mov_plate,self.fix_plate,min_age,max_age)


#############################################Rotation Object#####################################################

class Rot(object):

    def __init__(self,lat,lon,w,age_lb,age_ub,s1lat=0,s1lon=0,s1w=0,s1age_lb=0,s1age_ub=0):
        self.lat = float(lat)
        self.lon = float(lon)
#        self.azi = azi #add azimuth rots for oriented objects
        self.w = float(w)
        self.age_lb = float(age_lb)
        self.age_ub = float(age_ub)
        if self.age_lb==self.age_ub: self.wr = 0
        else: self.wr = w/(age_ub-age_lb)

        self.vlat = s1lat**2
        self.vlon = s1lon**2
        self.vw = s1w**2
        self.vage_lb = s1age_lb**2
        self.vage_ub = s1age_ub**2
        self.vwr = 0

        self.calc_matrix()

    #################Casting#################

    def to_matrix(self):
        try: return self.R
        except AttributeError: self.calc_matrix(); return self.R

    def to_list(self):
        return [self.lat,self.lon,self.w]

    def to_dict(self,add_data={}):
        rdict = {'lat':self.lat,'lon':self.lon,'rot':self.w,'rot_rate':self.wr,'start_age':self.age_lb, \
                 'stop_age':self.age_ub}
        return dict(rdict, **add_data)

    def copy(self):
        return Rot(self.lat,self.lon,self.w,self.age_lb,self.age_ub)

    def reverse_time(self):
        return Rot(self.lat,self.lon,-self.w,self.age_ub,self.age_lb)

    def backwards_rot(self):
        if self.age_lb>self.age_ub: return self.reverse_time()
        else: return self.copy()

    def forward_rot(self):
        if self.age_lb<self.age_ub: return self.reverse_time()
        else: return self.copy()

    #################Calculation#################

    def calc_matrix(self):
        rad_rot_angle = deg2rad(self.w)
        wx,wy,wz = latlon2cart(self.lat,self.lon)

        R1 = [wx*wx*(1-cos(rad_rot_angle))+cos(rad_rot_angle),
              wx*wy*(1-cos(rad_rot_angle))-wz*sin(rad_rot_angle),
              wx*wz*(1-cos(rad_rot_angle))+wy*sin(rad_rot_angle)]

        R2 = [wy*wx*(1-cos(rad_rot_angle))+wz*sin(rad_rot_angle),
              wy*wy*(1-cos(rad_rot_angle))+cos(rad_rot_angle),
              wy*wz*(1-cos(rad_rot_angle))-wx*sin(rad_rot_angle)]

        R3 = [wz*wx*(1-cos(rad_rot_angle))-wy*sin(rad_rot_angle),
              wz*wy*(1-cos(rad_rot_angle))+wx*sin(rad_rot_angle),
              wz*wz*(1-cos(rad_rot_angle))+cos(rad_rot_angle)]

        self.R = array([R1,R2,R3])

    def vel_point(self,lat,lon):
        arc_dis = Geodesic.WGS84.Inverse(lat,lon,self.lat,self.lon)['a12']
        if arc_dis>90: arc_dis = 180-arc_dis
        return self.wr*sin(deg2rad(arc_dis))

    def rotate(self,lat,lon,azi=0,d=1):

        geo_dict = Geodesic.WGS84.ArcDirect(lat,lon,azi,d)

        tmp_lat,tmp_lon = geo_dict['lat2'],geo_dict['lon2']

        nlat,nlon = self.rotate_site(lat,lon)
        ntmp_lat,ntmp_lon = self.rotate_site(tmp_lat,tmp_lon)

        geo_dict = Geodesic.WGS84.Inverse(nlat,nlon,ntmp_lat,ntmp_lon)

        return round(nlat,3),round(nlon,3),round(geo_dict['azi1'],3)

    def rotate_site(self,lat,lon):

        A = self.to_matrix()
        c = array(latlon2cart(lat,lon))
        nc = A @ c

        nlat,nlon = cart2latlon(*nc)

        return nlat,nlon

    #################Equality#################

    def eqrot(self,other_rot):
        rot_bool = self.w==other_rot.w and self.lat==other_rot.lat and self.lon==other_rot.lon and self.age_lb==other_rot.age_lb and self.age_ub==other_rot.age_ub
        cov_bool = self.vw==other_rot.vw and self.vlat==other_rot.vlat and self.vlon==other_rot.vlon and self.vage_lb==other_rot.vage_lb and self.vage_ub==other_rot.vage_ub
        return rot_bool and cov_bool

    #################Visualizing#################

    def plot_pole(self,m=None,**kwargs):
        return psk.plot_pole(self.lon,self.lat,0.0,sqrt(self.v1lon),sqrt(self.v1lat),m=m,**kwargs)

    #################Magic Functions#################

    ###########Arithmatic###########

    def __add__(self,other_rot):
        if not isinstance(other_rot, Rot): raise TypeError("Rotation addition is only defined with other Rot objects was given %s"%str(other_rot))
        R1 = self.to_matrix()
        R2 = other_rot.to_matrix()
        T = R2 @ R1
        return matrix2rot(T,self.age_lb,other_rot.age_ub)

    def __sub__(self,other_rot):
        return self+-other_rot

    def __neg__(self):
        return Rot(self.lat, self.lon, -self.w, self.age_lb, self.age_ub)

    def __invert__(self):
        if self.lon>0: nlon = self.lon-180
        else: nlon = self.lon+180
        return Rot(-self.lat,nlon,-self.w,self.age_lb,self.age_ub)

    ###########Equality###########

    def __hash__(self):
        return hash((self.age_lb,self.age_ub))

    def __eq__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return self.age_lb==other_rot.age_lb and self.age_ub==other_rot.age_ub

    def __ne__(self,other_rot):
        return (not self.__eq__(other_rot))

    def __lt__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return max(self.age_lb,self.age_ub)<max(other_rot.age_lb,other_rot.age_ub)

    def __gt__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return max(self.age_lb,self.age_ub)>max(other_rot.age_lb,other_rot.age_ub)

    def __le__(self,other_rot):
        return self.__lt__(other_rot) or self.__eq__(other_rot)

    def __ge__(self,other_rot):
        return self.__gt__(other_rot) or self.__eq__(other_rot)

    ###########Indexing & Sequencing###########

    def __getitem__(self,_slice):
        #check types and parse input
        if isinstance(_slice,slice): start,stop,step = _slice.start,_slice.stop,_slice.step
        elif isinstance(_slice,float) or isinstance(_slice,int): start,stop,step = _slice,None,None
        elif isinstance(_slice,list) or isinstance(_slice,tuple) or isinstance(_slice,ndarray):
            return [self.__getitem__(start) for start in _slice]
        else: raise TypeError("A slice of a rotation is not defined for %s"%str(_slice))

        #check boundedness of problem
        min_age,max_age = min(self.age_lb,self.age_ub),max(self.age_lb,self.age_ub)
        if start < min_age: raise ValueError("Start age for slice of roation less than rotation rotation was given %.3f but age bounds are %.3f-%.3f"%(float(start),self.age_lb,self.age_ub))
        if start > max_age: print("Indexed age not conatined in rotation extending roation by %.3f"%abs(start-max_age))
        if stop!=None and stop not in self: print("End age for slice of roation not conatined in rotation extending roation by %.3f"%abs(stop-max_age))

        if stop == None and step == None: #Index case (assumes you want the reconstruction pole)
            return Rot(self.lat,self.lon,self.wr*abs(self.age_lb-start),self.age_lb,start)
        elif stop != None and step == None: #Slice without step case (does not need to return reconstruction pole)
            return Rot(self.lat,self.lon,self.wr*(stop-start),start,stop)
        else: #Slice with step case
            rots = []
            for age in range(start,stop+step,step):
                if age==self.age_lb: continue
                rots.append(Rot(self.lat,self.lon,self.wr*abs(self.age_lb-age),self.age_lb,age))
            return rots

    def __contains__(self,age):
        try: age = float(age)
        except TypeError: raise TypeError("Expected type castable to float for age when checking if age in rotation got %s"%str(age))
        return age>=min(self.age_lb,self.age_ub) and age<=max(self.age_lb,self.age_ub)

    ###########Visualizing###########

    def __str__(self):
        return "lat = %.3f lon = %.3f w = %.3f\nage = %.3f-%.3f"%(self.lat,self.lon,self.w,float(self.age_lb),float(self.age_ub))

#############################################Custom Exceptions#####################################################

class RotationError(Exception):

    def __init__(self, message):
        super(RotationError, self).__init__(message)

#############################################Helper Functions#####################################################

def matrix2rot(R,age_lb,age_ub):
    lon = rad2deg(arctan2((R[0,2]-R[2,0]),(R[2,1]-R[1,2])))

    lat_num = R[1,0]-R[0,1]
    lat_den = sqrt((R[2,1]-R[1,2])**2 + (R[0,2]-R[2,0])**2 + (R[1,0]-R[0,1])**2)
    if lat_den==0: lat = 0
    else: lat = rad2deg(arcsin(lat_num/lat_den))

    w_num = lat_den
    w_den = R[0,0]+R[1,1]+R[2,2]-1
    w = rad2deg(arctan2(w_num,w_den))

    return Rot(lat,lon,w,age_lb,age_ub)

def latlon2cart(lat,lon):
    radlat = deg2rad(lat)
    radlon = deg2rad(lon)

    wx = cos(radlat)*cos(radlon)
    wy = cos(radlat)*sin(radlon)
    wz = sin(radlat)

    return [wx,wy,wz]

def cart2latlon(wx,wy,wz):
    lon = rad2deg(arctan2(wy,wx))
    lat = rad2deg(arctan2(wz,sqrt(wx**2 + wy**2)))

    return [lat,lon]

def merge_lists_without_duplicates(l1,l2):
    in_l1 = set(l1)
    in_l2 = set(l2)
    in_l2_not_in_l1 = in_l2 - in_l1
    return l1 + list(in_l2_not_in_l1)

def str_rots(rots):
    return reduce(lambda x,y: x+'\n'+y, map(str,rots))
