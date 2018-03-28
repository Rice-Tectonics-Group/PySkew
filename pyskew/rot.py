from numpy import *
from geographiclib.geodesic import Geodesic

class Rot(object):

    def __init__(self,lat,lon,w,age_lb,age_ub):
        self.lat = lat
        self.lon = lon
#        self.azi = azi #add azimuth rotations for oriented objects
        self.w = w
        self.age_lb = age_lb
        self.age_ub = age_ub
        self.wr = w/(age_ub-age_lb)
        self.calc_matrix()

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

    def to_matrix(self):
        try: return self.R
        except AttributeError: self.calc_matrix(); return self.R

    def to_list(self):
        return [self.lat,self.lon,self.w]

    def __add__(self,other_rot):
        if not isinstance(other_rot, Rot): raise TypeError("Rotation addition is only defined with other Rot objects please pass in a Rot object to preform addition with this object")
        R1 = self.to_matrix()
        R2 = other_rot.to_matrix()
        T = R2.dot(R1)
        return matrix2rot(T,self.age_lb,other_rot.age_ub)

    def __neg__(self):
        return Rot(self.lat, self.lon, -self.w, self.age_lb, self.age_ub)

    def __sub__(self,other_rot):
        return self+-other_rot

    def plot_pole(self,m,**kwargs):
        m.scatter(self.lon,self.lat,**kwargs)

    def antipode(self):
        if self.lon>0: nlon = self.lon-180
        else: nlon = self.lon+180
        return Rot(-self.lat,nlon,-self.w,self.age_lb,self.age_ub)

    def reverse_time(self):
        return Rot(self.lat,self.lon,-self.w,self.age_ub,self.age_lb)

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
        nc = A.dot(c)

        nlat,nlon = cart2latlon(*nc)

        return nlat,nlon

    def __str__(self):
        return "lat=%.3f lon=%.3f w=%.3f\nage = %.3f-%.3f"%(self.lat,self.lon,self.w,float(self.age_lb),float(self.age_ub))

def matrix2rot(R,age_lb,age_ub):
    lon = rad2deg(arctan2((R[0,2]-R[2,0]),(R[2,1]-R[1,2])))

    lat_num = R[1,0]-R[0,1]
    lat_den = sqrt((R[2,1]-R[1,2])**2 + (R[0,2]-R[2,0])**2 + (R[1,0]-R[0,1])**2)
    lat = rad2deg(arcsin(lat_num/lat_den))

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
