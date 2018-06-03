from numpy import *
from numpy.linalg import svd
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

    def __init__(self,mov_plate,fix_plate,rots,comment=None,sources={},parent=None):
        self.mov_plate = str(mov_plate)
        self.fix_plate = str(fix_plate)
        if comment==None: comment = "Motion of %s relative to %s"%(mov_plate,fix_plate)
        self.comment = str(comment)
        self.sources = sources
        if self.sources=={}: self.sources = {rot:"" for rot in rots}
        if any(list(map(lambda x: not isinstance(x,Rot), rots))):
            raise ValueError("All values in rotsn argument must be Rot objects, was given %s"%str(rots))
        try: tmp_rots = [rot.backwards_rot() for rot in rots]
        except ValueError: raise ValueError("rots argument must be list-like, was given %s"%str(rots))
        tmp_rots.sort() #insure self.rots is always sorted by age

        self.parent = parent
        self.children = []
        if not isinstance(parent,type(None)): self.rot,self.age = tmp_rots[0],tmp_rots[0].age_f
        else: self.rot,self.age = None,tmp_rots[0].age_i
        for i,rot in enumerate(tmp_rots):
            if rot.age_i!=self.age: continue
            chrots = [rot] + [chrot for chrot in tmp_rots if chrot.age_i==rot.age_f]
            if i+1<len(tmp_rots) and (any([chrot.age_f>tmp_rots[(i+1)].age_f for chrot in chrots]) or any([chrot==tmp_rots[(i+1)] for chrot in chrots])):
                raise ValueError("One of your stage rotation is providing duplicate information to your reconstruction rotation\n%s"%(str(tmp_rots[i+1])))
            chsources = {self.sources[chrot]:chrot for chrot in chrots if chrot in self.sources}
            self.children.append(PlateReconstruction(self.mov_plate,self.fix_plate,chrots,comment=self.comment,sources=chsources,parent=self))

    #################I/O#################

    @classmethod
    def read_csv(self,fin,sep='\t'):
        with open(fin,'r') as f:
            mov_plate,fix_plate,comment = f.readline().rstrip('\n').split('\t')[:3]
        df = pd.read_csv(fin,header=1,sep='\t')
        df.fillna(0, inplace=True)
        errors = [{'cov':vs_to_cov(*latlonrot2cart(lat,lon,w,zeros([3,3]))[0],vs),'s1age_f':s1age_f,'s1age_i':s1age_i} for lat,lon,w,vs,s1age_i,s1age_f in zip(df['lat'],df['lon'],df['rot'],df[['vlat','vlon','vrot','vlatlon','vlatrot','vlonrot']].values.tolist(),df['s1start_age'],df['s1stop_age'])]
        rots_list = df[['lat','lon','rot','start_age','stop_age']].values.tolist()
        rots = [Rot(*rot,**err) for rot,err in zip(rots_list,errors)]
        sources = df['source'].tolist()
        new_sources = {rot:source for source,rot in zip(sources,rots)}
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
        df = pd.DataFrame([(r.backwards_rot()).to_dict(add_data={'source':(self.sources[r] if r in self.sources.keys() else "")}) for r in self.get_rots()])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','s1start_age','s1stop_age','vlat','vlon','vrot','vlatlon','vlatrot','vlonrot','source']]
        return df

    def reconst_to_df(self):
        df = pd.DataFrame([(r).to_dict(add_data={'source':'this study'}) for r in self.get_rrots()])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','s1start_age','s1stop_age','vlat','vlon','vrot','vlatlon','vlatrot','vlonrot','source']]
        return df

    def stages_to_df(self):
        df = pd.DataFrame([(r).to_dict(add_data={'source':'this study'}) for r in self.get_srots()])
        df = df[['lat','lon','rot','rot_rate','start_age','stop_age','s1start_age','s1stop_age','vlat','vlon','vrot','vlatlon','vlatrot','vlonrot','source']]
        return df

    def list_ages(self):
        """
        Return ages for which there are rotations to the age of the root
        """
        ages=[]
        for rot in self.get_rrots():
            if rot.age_i not in ages: ages.append(rot.age_i)
            if rot.age_f not in ages: ages.append(rot.age_f)
        if self.age in ages: ages.remove(self.age)
        return ages

    def get_age_bounds(self):
        ages = [self.age]+self.list_ages()
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

    def get_rots(self):
        if self.rot==None: rots = []
        else: rots = [self.rot]
        for child in self.children:
            for rot in child.get_rots():
                rots.append(rot)
        return rots

    def get_srots(self):
        if self.rot==None: rots = []
        else: rots = [self.rot]
        if len(self.children)==0: return rots
        rots += self.children[0].get_srots()
        for child in self.children[1:]:
            for rot in child.get_srots():
                rots.append(self[rots[-1].age_f:self.age]+rot) #THIS IS WRONG ADDS STAGE POLE TO RECONSTRUCTION POLES
        return rots

    def get_rrots(self):
        if self.rot==None:
            rots = []
            for child in self.children:
                for rot in child.get_rrots():
                    rots.append(rot)
        else:
            rots = [self.rot]
            for child in self.children:
                for rot in child.get_rrots():
                    rots.append(self.rot+rot)
        return rots

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

    def complex__add__(self,other_reconst):
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

    def __add__(self,other_reconst):
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
        new_rots = [-r for r in self.get_rots()]
        return PlateReconstruction(self.fix_plate,self.mov_plate,new_rots,comment=self.comment,sources=self.sources)

    def __invert__(self):
        new_rots = [~r for r in self.get_rots()]
        return PlateReconstruction(self.mov_plate,self.fix_plate,new_rots,comment=self.comment,sources=self.sources)

    ###########Indexing & Sequencing###########

    def __len__(self):
        return len(self.children)
######################################################THIS SECTION NEEDS OVERHAUL##################################
    def __getitem__(self,_slice):
        #check types and parse input
        if isinstance(_slice,slice): start,stop,step = _slice.start,_slice.stop,_slice.step
        elif isnum(_slice): start,stop,step = float(_slice),None,None
        elif isiter(_slice):
            return [self.__getitem__(slice(*item)) if isiter(item) else self.__getitem__(item) for item in _slice]
        else: raise TypeError("A slice of a rotation is not defined for %s"%str(_slice))

        #check boundedness of problem
        min_age,max_age = self.get_age_bounds()
        if start < min_age:
            raise ValueError("Start age of %.3f is less than youngest reconstruction for %s relative to %s"%(float(start),self.mov_plate,self.fix_plate))
        if start > max_age:
            print("Start age %.3f greater than oldest reconstruction for %s relative to %s, extending final rotation"%(float(start),self.mov_plate,self.fix_plate))
            return self.get_rrots()[-1][start]

        start_parent,prev=self,None
        for child in self.children:
            if start==child.age: start_rot = child.rot; break
            elif child.age>start:
                if prev==None: start_rot = child.rot[start]
                else: start_rot = prev.rot+(prev.rot.reverse_time()+child.rot)[start]; break
            elif child.get_age_bounds()[1]>start: start_rot,start_parent = child[start],child; break #will be reconsturction to child
            else: prev=child

        if stop==None:
            if start==self.age: return Rot()
            elif start_parent==self: return start_rot
            else: return start_parent.rot+start_rot
        elif step==None:
            return start_rot.forward_rot() + start_parent[stop]
        else:
            rots = []
            for next in arange(start+step,stop+step,step):
                rots.append(self[start:next])
            return rots

        raise RuntimeError("Couldn't find a rotation or a problem with input in __getitem__ for reconstruction %s this is almost certainly a bug, contact a developer at the github page."%(str(self)))

    def __setitem__(self,idx,value):
        if isinstance(value,Rot): nrot=value
        elif isinstance(idx,slice) and len(value)==3: nrot = Rot(*value,idx.start,idx.stop)
        else: raise TypeError("Cannot add roation to PlateReconstruction for index: %s and value: %s. You must provide either a Rot object for value and any index (index will be ignored) or a slice for index (i.e. [0:10]) which represents your age bounds and a list with values [lat,lon,w] for your value."%(str(idx),str(value)))
        if nrot.age_i==self.age:
            for child in self.children:
                if nrot==child.rot: self.children.remove(child); del self.sources[child.rot]
            self.children.append(PlateReconstruction(self.mov_plate,self.fix_plate,[nrot],comment=self.comment,parent=self))
            self.sources[nrot] = ""
        elif nrot.age_i>self.get_age_bounds()[1] and self.children!=[]:
            oldest = self.children[-1]
            if oldest.children!=[]: oldest.__setitem__(idx,nrot)
            else:
                nrot_interp = oldest.rot[nrot.age_i]+nrot
                self.children.append(PlateReconstruction(self.mov_plate,self.fix_plate,[nrot_interp],comment=self.comment,parent=self))
                self.sources[nrot_interp] = ""
        else:
            for child in self.children:
                child.__setitem__(idx,nrot)

    def __delitem__(self,value):
        if isinstance(value,Rot):
            for child in self.children:
                if child.rot==Rot: self.children.remove(child); del self.sources[child.rot]
                elif value in child: child.__delitem__(value)
            else: raise ValueError("There is no equivlent rotation in this reconstruction cannot remove:\n%s"%str(value))
        elif isnum(value):
            for child in self.children:
                if value in child.rot: self.children.remove(child); del self.sources[child.rot]
                elif value in child: child.__delitem__(value)
        else:
            raise TypeError("There is no way to remove %s from reconstruction"%str(value))

    def __contains__(self,value):
        if isnum(value):
            return any([(value in rot) for rot in self.get_rots()])
        elif isinstance(value,Rot):
            return ((value in self.get_rots()))
        else:
            raise TypeError("cannot check if %s is in a reconstruction object"%str(value))

    ###########Visualizing###########

    def __str__(self):
        min_age,max_age = self.get_age_bounds()
        return "Reconstruction Rotations for %s relative to %s, from %.2f to %.2f"%(self.mov_plate,self.fix_plate,min_age,max_age)

    def print_tree(self):
        print(self.tree_str())

    def tree_str(self):
        out_str = str(self.age)
        if self.children: out_str += "\n|\nv\n"
        for child in self.children:
            if child!=self.children[-1]: out_str += child.tree_str() + ","
            else: out_str += child.tree_str()
        return out_str


#############################################Rotation Object#####################################################

class Rot(object):

    def __init__(self,lat=0,lon=0,w=0,age_i=0,age_f=0,s1lat=0,s1lon=0,s1w=0,s1age_i=0,s1age_f=0,cov=None):
        self.lat = float(lat)
        self.lon = float(lon)
#        self.azi = azi #add azimuth rots for oriented objects
        self.w = float(w)
        self.age_i = float(age_i)
        self.age_f = float(age_f)
        if self.age_i==self.age_f: self.wr = 0
        else: self.wr = w/(age_f-age_i)

        if isinstance(cov,type(None)):
            self.cov = diag([s1lat**2,s1lon**2,s1w**2])
        else: self.cov = cov
        self.vage_i = s1age_i**2
        self.vage_f = s1age_f**2
        if self.w==0 or (self.age_i==0 and self.age_f==0): self.vwr=0
        else: self.vwr = (self.wr**2)*((self.cov[2,2]/(self.w**2)) + ((self.vage_f+self.vage_i)/(self.age_f+self.age_i)**2)) #IGNORES COVARIANCE HERE

        self.calc_matrix()

    #################Casting#################

    def to_matrix(self):
        try: return self.R
        except AttributeError: self.calc_matrix(); return self.R

    def get_matrix_cov(self):
        try: return self.R_cov
        except AttributeError: self.calc_matrix(); return self.R_cov

    def to_list(self):
        return [self.lat,self.lon,self.w]

    def to_dict(self,add_data={}):
        vs = cov_to_vs(self.lat,self.lon,self.w,self.cov)
        rdict = {'lat':self.lat,'lon':self.lon,'rot':self.w,'rot_rate':self.wr,'start_age':self.age_i, \
                 'stop_age':self.age_f,'s1start_age':sqrt(self.vage_i),'s1stop_age':sqrt(self.vage_f),
                 'vlat':vs[0],'vlon':vs[1],'vrot':vs[2],'vlatlon':vs[3],'vlatrot':vs[4],'vlonrot':vs[5]}
        return dict(rdict, **add_data)

    def copy(self):
        return Rot(self.lat,self.lon,self.w,self.age_i,self.age_f,s1age_f=sqrt(self.vage_f),s1age_i=sqrt(self.vage_i),cov=self.cov)

    def reverse_time(self):
#        ncov = self.inv_wcovs(self.cov)
        return Rot(self.lat,self.lon,-self.w,self.age_f,self.age_i,s1age_f=sqrt(self.vage_i),s1age_i=sqrt(self.vage_f),cov=self.cov)

    def backwards_rot(self):
        if self.age_i>self.age_f: return self.reverse_time()
        else: return self.copy()

    def forward_rot(self):
        if self.age_i<self.age_f: return self.reverse_time()
        else: return self.copy()

    def isnull(self): #any matrix which does not rotate is considered null
        return self.w==0

    #################Calculation#################

    def calc_matrix(self):
        rad_rot_angle = deg2rad(self.w)
        (wx,wy,wz),cart_cov = latlon2cart(self.lat,self.lon,self.cov) #also propogates errors for rotation angle

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

        J = array([[2*wx*(1-cos(rad_rot_angle)),0,0,wx*wx*sin(rad_rot_angle)-sin(rad_rot_angle)],
                    [wy*(1-cos(rad_rot_angle)),wx*(1-cos(rad_rot_angle)),-sin(rad_rot_angle),wx*wy*sin(rad_rot_angle)-wz*cos(rad_rot_angle)],
                    [wz*(1-cos(rad_rot_angle)),sin(rad_rot_angle),wx*(1-cos(rad_rot_angle)),wx*wz*sin(rad_rot_angle)+wy*cos(rad_rot_angle)],
                    [wy*(1-cos(rad_rot_angle)),wx*(1-cos(rad_rot_angle)),sin(rad_rot_angle),wy*wx*sin(rad_rot_angle)+wz*cos(rad_rot_angle)],
                    [0,2*wy*(1-cos(rad_rot_angle)),0,wy*wy*sin(rad_rot_angle)-sin(rad_rot_angle)],
                    [-sin(rad_rot_angle),wz*(1-cos(rad_rot_angle)),wy*(1-cos(rad_rot_angle)),wy*wz*sin(rad_rot_angle)-wx*cos(rad_rot_angle)],
                    [wz*(1-cos(rad_rot_angle)),-sin(rad_rot_angle),wx*(1-cos(rad_rot_angle)),wz*wx*sin(rad_rot_angle)-wy*cos(rad_rot_angle)],
                    [sin(rad_rot_angle),wz*(1-cos(rad_rot_angle)),wy*(1-cos(rad_rot_angle)),wz*wy*sin(rad_rot_angle)+wx*cos(rad_rot_angle)],
                    [0,0,2*wz*(1-cos(rad_rot_angle)),wz*wz*sin(rad_rot_angle)-sin(rad_rot_angle)]])

        self.R_cov = J @ cart_cov @ J.T

    def vel_point(self,lat,lon):
        arc_dis = Geodesic.WGS84.Inverse(lat,lon,self.lat,self.lon)['a12'] #SHOULD HAVE ARCDIS ERROR AS LENGTH OF ELLIPSE AROUND POLE IN AZI DIRECTION
        if arc_dis>90: arc_dis = 180-arc_dis
        return self.wr*sin(deg2rad(arc_dis))

    #NEED TO CHANGE THE METHOD HERE AZI AND PHI REDUNDANT NOW ERROR PROP WORKS WE DON'T NEED AZI
    def rotate(self,lat,lon,azi=0,a=None,b=None,phi=None,cov=None,d=1):

        if isinstance(cov,type(None)):
            if isinstance(a,type(None)) or isinstance(b,type(None)) or isinstance(phi,type(None)):
                cov = zeros([2,2])
            else:
                cov = ellipse_to_cov(a,b,phi)

        geo_dict = Geodesic.WGS84.ArcDirect(lat,lon,azi,d)

        tmp_lat,tmp_lon = geo_dict['lat2'],geo_dict['lon2']

        nlat,nlon,ell,ncov = self.rotate_site(lat,lon,cov=cov)
        ntmp_lat,ntmp_lon,tmp_ell,ntmp_cov = self.rotate_site(tmp_lat,tmp_lon,cov=cov)

        geo_dict = Geodesic.WGS84.Inverse(nlat,nlon,ntmp_lat,ntmp_lon)

        return round(nlat,3),round(nlon,3),round(geo_dict['azi1'],3),ell #I think this already takes into account the azimuth uncertainty but I'm not sure

    def rotate_site(self,lat,lon,cov=None):

        if isinstance(cov,type(None)):
            cov = zeros([2,2])

        A = self.to_matrix()
        cov_A = self.get_matrix_cov()
        c,cart_cov = latlon2cart(lat,lon,cov)
        nc = A @ c
        Jr = array([[c[0],0,0,c[1],0,0,c[2],0,0],[0,c[0],0,0,c[1],0,0,c[2],0],[0,0,c[0],0,0,c[1],0,0,c[2]]])
        ncart_cov = A @ cart_cov @ A.T + Jr @ cov_A @ Jr.T

        (nlat,nlon),ncov = cart2latlon(*nc,ncart_cov)

        return nlat,nlon,cov_to_ellipse(ncov),ncov

    #################Equality#################

    def eqrot(self,other_rot):
        rot_bool = self.w==other_rot.w and self.lat==other_rot.lat and self.lon==other_rot.lon and self.age_i==other_rot.age_i and self.age_f==other_rot.age_f
        return rot_bool

    #################Visualizing#################

    def plot_pole(self,m=None,**kwargs):
        return psk.plot_pole(self.lon,self.lat,0.0,sqrt(self.v1lon),sqrt(self.v1lat),m=m,**kwargs)

    #################Magic Functions#################

    ###########Arithmatic###########

    @staticmethod
    def get_R1J(R2):
        return hstack([vstack([R2.T,zeros([6,3])]),vstack([zeros([3,3]),R2.T,zeros([3,3])]),vstack([zeros([6,3]),R2.T])])

    @staticmethod
    def get_R2J(R1):
        rows = []
        for row in R1.T:
            col = []
            for r in row:
                col.append(eye(3)*r)
            rows.append(col)
        return hstack([vstack(col) for col in rows])

    @staticmethod
    def old_get_R1J(R2):
        Js = []
        for r in R2.T:
            Js.append(hstack([vstack([r,zeros([2,3])]),vstack([zeros([1,3]),r,zeros([1,3])]),vstack([zeros([2,3]),r])]))
        return vstack(Js)

    @staticmethod
    def old_get_R2J(R1):
        Js = []
        for r in R1.T:
            r = array([r])
            Js.append(vstack([hstack([r.T,zeros([3,2])]),hstack([zeros([3,1]),r.T,zeros([3,1])]),hstack([zeros([3,2]),r.T])]))
        return hstack(Js)

    def __add__(self,other_rot):
        if not isinstance(other_rot, Rot): raise TypeError("Rotation addition is only defined with other Rot objects was given %s"%str(other_rot))
        R1 = self.to_matrix()
        R1_cov = self.get_matrix_cov()
        R2 = other_rot.to_matrix()
        R2_cov = other_rot.get_matrix_cov()
        T = R2 @ R1
        J1 = self.get_R1J(R2.T)
        J2 = self.get_R2J(R1.T)
        cov_T = J1 @ R1_cov @ J1.T + J2 @ R2_cov @ J2.T

        return matrix2rot(T,self.age_i,other_rot.age_f,cov=cov_T)

    @staticmethod
    def inv_wcovs(cov):
        ncov = -cov.copy()
        ncov[0,0],ncov[1,1],ncov[2,2] = -ncov[0,0],-ncov[1,1],-ncov[2,2]
        return ncov

    def __sub__(self,other_rot):
        return self+-other_rot

    def __neg__(self):
#        ncov = self.inv_wcovs(self.cov)
        return Rot(self.lat, self.lon, -self.w, self.age_i, self.age_f, s1age_f=sqrt(self.vage_f), s1age_i=sqrt(self.vage_i), cov=self.cov)

    def __invert__(self):
        if self.lon>0: nlon = self.lon-180
        else: nlon = self.lon+180
        return Rot(-self.lat,nlon,-self.w,self.age_i,self.age_f, s1age_f=sqrt(self.vage_f), s1age_i=sqrt(self.vage_i), cov=self.cov)

    ###########Equality###########

    def __hash__(self):
        return hash((self.age_i,self.age_f))

    def __eq__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return self.age_i==other_rot.age_i and self.age_f==other_rot.age_f

    def __ne__(self,other_rot):
        return (not self.__eq__(other_rot))

    def __lt__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return max(self.age_i,self.age_f)<max(other_rot.age_i,other_rot.age_f)

    def __gt__(self,other_rot):
        if not isinstance(other_rot, Rot): return False
        return max(self.age_i,self.age_f)>max(other_rot.age_i,other_rot.age_f)

    def __le__(self,other_rot):
        return self.__lt__(other_rot) or self.__eq__(other_rot)

    def __ge__(self,other_rot):
        return self.__gt__(other_rot) or self.__eq__(other_rot)

    ###########Indexing & Sequencing###########

    def __getitem__(self,_slice): #MAY NEED TO CONSIDER INTERPOLATION'S CONTRIBUTIN TO ERROR
        #check types and parse input
        if isinstance(_slice,slice): start,stop,step = _slice.start,_slice.stop,_slice.step
        elif isnum(_slice): start,stop,step = float(_slice),None,None
        elif isiter(_slice):
            return [self.__getitem__(start) for start in _slice]
        else: raise TypeError("A slice of a rotation is not defined for %s"%str(_slice))

        #check boundedness of problem
        min_age,max_age = min(self.age_i,self.age_f),max(self.age_i,self.age_f)
        if start < min_age: raise ValueError("Start age for slice of roation less than rotation rotation was given %.3f but age bounds are %.3f-%.3f"%(float(start),self.age_i,self.age_f))
        if start > max_age: print("Indexed age not conatined in rotation extending roation by %.3f"%abs(start-max_age))
        if stop!=None and stop not in self: print("End age for slice of roation not conatined in rotation extending roation by %.3f"%abs(stop-max_age))

        if stop == None and step == None: #Index case (assumes you want the reconstruction pole)
            cov = array(self.cov)
            if self.wr!=0 and self.age_i!=0: cov[2,2] = (self.w**2)*((self.vwr/(self.wr**2)) + (self.vage_i/(self.age_i**2)))
            elif self.wr!=0 and self.age_i==0: cov[2,2] = (self.w**2)*(self.vwr/(self.wr**2))
            return Rot(self.lat,self.lon,self.wr*abs(self.age_i-start),self.age_i,start,cov=cov)
        elif stop != None and step == None: #Slice without step case (does not need to return reconstruction pole)
            cov = array(self.cov)
            if self.wr!=0: cov[2,2] = (self.w**2)*((self.vwr/(self.wr**2)))
            return Rot(self.lat,self.lon,self.wr*(stop-start),start,stop,cov=cov)
        else: #Slice with step case
            rots = []
            for age in range(start,stop+step,step):
                if age==self.age_i: continue
                cov = array(self.cov)
                if self.wr!=0 and self.age_i!=0: cov[2,2] = (self.w**2)*((self.vwr/(self.wr**2)) + (self.vage_i/(self.age_i**2)))
                elif self.wr!=0 and self.age_i==0: cov[2,2] = (self.w**2)*(self.vwr/(self.wr**2))
                rots.append(Rot(self.lat,self.lon,self.wr*abs(self.age_i-age),self.age_i,age,cov=cov))
            return rots

    def __contains__(self,age):
        try: age = float(age)
        except TypeError: raise TypeError("Expected type castable to float for age when checking if age in rotation got %s"%str(age))
        return age>=min(self.age_i,self.age_f) and age<=max(self.age_i,self.age_f)

    ###########Visualizing###########

    def __str__(self):
        return "lat = %.3f lon = %.3f w = %.3f\nage = %.3f-%.3f"%(self.lat,self.lon,self.w,float(self.age_i),float(self.age_f))

#############################################Custom Exceptions#####################################################

class RotationError(Exception):

    def __init__(self, message):
        super(RotationError, self).__init__(message)

#############################################Helper Functions#####################################################

def matrix2rot(R,age_i,age_f,cov=None):
    #May not need this consider removing later
#    R = array([R[i,j] if R[i,j]>1e-10 else 0 for i in range(R.shape[0]) for j in range(R.shape[1])]).reshape(R.shape)

    lon = rad2deg(arctan2((R[0,2]-R[2,0]),(R[2,1]-R[1,2])))

    lat_num = R[1,0]-R[0,1]
    lat_den = sqrt((R[2,1]-R[1,2])**2 + (R[0,2]-R[2,0])**2 + (R[1,0]-R[0,1])**2)
    if lat_den==0: lat = 0
    else: lat = rad2deg(arcsin(lat_num/lat_den))

    w_num = lat_den
    w_den = R[0,0]+R[1,1]+R[2,2]-1
    w = rad2deg(arctan2(w_num,w_den))

    if R.all()==0 or (R==eye(3)).all(): J=eye(3,9)
    else:
        lon_arg = (R[0,2]-R[2,0])/(R[2,1]-R[1,2])
        der_asin = lambda x: 1/sqrt(1-x**2)
        der_atan = lambda x: 1/(1+x**2)
        J = array([[0,der_asin(lat_num/lat_den)*(((R[2,0]-R[0,2])**2 + (R[2,1]-R[1,2])**2)/((R[1,0]-R[0,1])**2 + (R[2,0]-R[0,2])**2 + (R[2,1]-R[1,2])**2)**(3/2)),der_asin(lat_num/lat_den)*(((R[1,0]-R[0,1])*(R[0,2]-R[2,0]))/((R[1,0]-R[0,1])**2 + (R[0,2]-R[2,0])**2 + (R[2,1]-R[1,2])**2)**(3/2)),
                    der_asin(lat_num/lat_den)*((-(R[2,0]-R[0,2])**2 - (R[2,1]-R[1,2])**2)/((R[1,0]-R[0,1])**2 + (R[2,0]-R[0,2])**2 + (R[2,1]-R[1,2])**2)**(3/2)),0,der_asin(lat_num/lat_den)*((R[1,0]-R[0,1])*(R[1,2]-R[2,1]))/((R[1,0]-R[0,1])**2 + (R[0,2]-R[2,0])**2 + (R[2,1]-R[1,2])**2)**(3/2),
                    der_asin(lat_num/lat_den)*(((R[1,0]-R[0,1])*(R[2,0]-R[0,2]))/((R[1,0]-R[0,1])**2 + (R[2,0]-R[0,2])**2 + (R[2,1]-R[1,2])**2)**(3/2)),der_asin(lat_num/lat_den)*(((R[1,0]-R[0,1])*(R[2,1]-R[1,2]))/((R[1,0]-R[0,1])**2 + (R[0,2]-R[2,0])**2 + (R[2,1]-R[1,2])**2)**(3/2)),0],

                   [0,0,-der_atan(lon_arg)/(R[2,1]-R[1,2]),
                    0,0,-der_atan(lon_arg)*((R[0,2]-R[2,0])/((R[2,1]-R[1,2])**2)),
                    der_atan(lon_arg)/(R[2,1]-R[1,2]),der_atan(lon_arg)*((R[0,2]-R[2,0])/((R[2,1]-R[1,2])**2)),0],

                   [der_atan(w_num/w_den)*(w_num/(w_den**2)),der_atan(w_num/w_den)*((R[1,0]-R[0,1])/(w_den*w_num)),-der_atan(w_num/w_den)*((R[0,2]-R[2,0])/(w_den*w_num)),
                    -der_atan(w_num/w_den)*((R[1,0]-R[0,1])/(w_den*w_num)),der_atan(w_num/w_den)*(w_num/(w_den**2)),der_atan(w_num/w_den)*((R[2,1]-R[1,2])/(w_den*w_num)),
                    der_atan(w_num/w_den)*((R[0,2]-R[2,0])/(w_den*w_num)),-der_atan(w_num/w_den)*((R[2,1]-R[1,2])/(w_den*w_num)),der_atan(w_num/w_den)*(w_num/(w_den**2))]])

    new_cov = (((180/pi)**2)*(J @ cov @ J.T))

    if isnan(new_cov).any(): print("ERROR IN COVARIANCE CALCULATION CAUSED BY THAT DAMN NULL ROTATION, CALL DEV")

    return Rot(lat,lon,w,age_i,age_f,cov=new_cov)

def latlonrot2cart(lat,lon,w,cov):
    radlat = deg2rad(lat)
    radlon = deg2rad(lon)
    radw = deg2rad(w)

    wx = radw*cos(radlat)*cos(radlon)
    wy = radw*cos(radlat)*sin(radlon)
    wz = radw*sin(radlat)

    J = array([[-radw*sin(radlat)*cos(radlon),-radw*sin(radlon)*cos(radlat),cos(radlat)*cos(radlon)],
               [-radw*sin(radlat)*sin(radlon),radw*cos(radlon)*cos(radlat),cos(radlat)*sin(radlon)],
               [radw*cos(radlat),0,sin(radlat)]])

    cart_cov = J @ (((pi/180)**2)*cov) @ J.T

    return [wx,wy,wz],cart_cov

def cart2latlonrot(wx,wy,wz,cart_cov):
    lat = rad2deg(arctan2(wz,sqrt(wx**2 + wy**2)))
    lon = rad2deg(arctan2(wy,wx))
    w = rad2deg(sqrt(wx**2 + wy**2 + wz**2))

    J = array([[(-wx*wz)/(sqrt(wx**2+wy**2)*(wx**2+wy**2+wz**2)),
               (-wy*wz)/(sqrt(wx**2+wy**2)*(wx**2+wy**2+wz**2)),
               sqrt(wx**2 + wy**2)/(wx**2 + wy**2 + wz**2)],
               [-wy/(wx**2+wy**2),(wx)/(wx**2+wy**2),0],
               [wx/deg2rad(w),wy/deg2rad(w),wz/deg2rad(w)]])

    latlonrot_cov = (((180/pi)**2)*(J @ cart_cov @ J.T))

    return [lat,lon,w],latlonrot_cov

def latlon2cart(lat,lon,cov):
    radlat = deg2rad(lat)
    radlon = deg2rad(lon)

    wx = cos(radlat)*cos(radlon)
    wy = cos(radlat)*sin(radlon)
    wz = sin(radlat)

    if len(cov)==3:
        J = array([[-sin(radlat)*cos(radlon),-sin(radlon)*cos(radlat),0],
                   [-sin(radlat)*sin(radlon),cos(radlon)*cos(radlat),0],
                   [cos(radlat),0,0],
                   [0,0,1]])
    else:
        J = array([[-sin(radlat)*cos(radlon),-sin(radlon)*cos(radlat)],
                   [-sin(radlat)*sin(radlon),cos(radlon)*cos(radlat)],
                   [cos(radlat),0]])

    cart_cov = J @ (((pi/180)**2)*cov) @ J.T

    return [wx,wy,wz],cart_cov

def cart2latlon(wx,wy,wz,cart_cov):
    lon = rad2deg(arctan2(wy,wx))
    lat = rad2deg(arctan2(wz,sqrt(wx**2 + wy**2)))

    J = array([[(-wx*wz)/(sqrt(wx**2+wy**2)*(wx**2+wy**2+wz**2)),
               (-wy*wz)/(sqrt(wx**2+wy**2)*(wx**2+wy**2+wz**2)),
               sqrt(wx**2 + wy**2)/(wx**2 + wy**2 + wz**2)],
               [-wy/(wx**2+wy**2),(wx)/(wx**2+wy**2),0]])

    latlon_cov = (((180/pi)**2)*(J @ cart_cov @ J.T))

    return [lat,lon],latlon_cov

def ellipse_to_cov(sa,sb,phi):
    rphi = deg2rad(phi)
    var1 = sa**2 * cos(rphi)**2 + sb**2 * sin(rphi)**2 #vlat
    var2 = sa**2 * sin(rphi)**2 + sb**2 * cos(rphi)**2 #vlon
    covar12 = (sa**2 - sb**2) * sin(rphi) * cos(rphi)
    return array([[var1,covar12],[covar12,var2]])

def cov_to_ellipse(cov): # NEED TO ACCOUNT FOR THE FACT THAT LAT AND LON HAVE DIFFERENT METRICS BAISED ON POSITION
    w,v = linalg.eig(cov)
    vmax = v[:,list(w).index(max(w))]
    phi = rad2deg(arctan2(vmax[1],vmax[0]))
    return tuple((*sorted(sqrt(w),reverse=True),phi))

def vs_to_cov(wx,wy,wz,vs): #I NEED TO CHANGE THE DATA READ WRITE FORMAT RIGHT NOW IT IS KOI SPECIFIC
    return cart2latlonrot(wx,wy,wz,array([[vs[0],vs[3],vs[4]],
                     [vs[3],vs[1],vs[5]],
                     [vs[4],vs[5],vs[2]]]))[-1]

def cov_to_vs(lat,lon,w,cov):
    _,cov = latlonrot2cart(lat,lon,w,cov)
    return [cov[0,0],cov[1,1],cov[2,2],cov[0,1],cov[0,2],cov[1,2]]

def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

def get_R_star(R,p=1):
    """
    Returns matrix of all multiplicative combinations of the elements of R raised to a power p
    """
    R_star = zeros([R.shape[0]*R.shape[1],R.shape[0]*R.shape[1]])
    for i,r1 in enumerate(R.flatten()):
        for j,r2 in enumerate(R.flatten()):
            if p<0 and r1*r2==0:
                R_star[i,j] = 0
            else:
                R_star[i,j] = ((r1*r2)**p)
    return R_star

def isnum(obj):
    try:
        float(obj)
        return True
    except (ValueError,TypeError) as e:
        return False

def isiter(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False

def merge_lists_without_duplicates(l1,l2):
    in_l1 = set(l1)
    in_l2 = set(l2)
    in_l2_not_in_l1 = in_l2 - in_l1
    return l1 + list(in_l2_not_in_l1)

def str_rots(rots):
    return reduce(lambda x,y: x+'\n'+y, map(str,rots))
