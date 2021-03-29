import os
import numpy as np
import rdp
from pyrot.rot import cart2latlon,latlon2cart
from geographiclib.geodesic import Geodesic

def find_corners(lat_lon_file,tolerance=1,min_angle=np.pi*0.22,geoid=Geodesic(6371.,0.),lsz=None):
    """
    Credit to unutbu from Stack Overflow
    (https://stackoverflow.com/questions/14631776/calculate-turning-points-pivot-points-in-trajectory-path)
    Modifed by PySkew author Kevin Gaastra to run on a sphere

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

    # Use the Ramer-Douglas-Peucker algorithm to simplify the path
    # http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
    # Python implementation: https://github.com/sebleier/RDP/
    simplified = np.array(rdp.rdp(points.tolist(), tolerance, dist=create_geo_gc_dist(geoid)))
#    simplified = np.array(kmg9_rdp(points.tolist(), tolerance, create_geo_gc_dist(geoid),lsz=lsz))

    print("simplified number of points = %d"%len(simplified))

    # compute the direction vectors on the simplified curve
#    directions = np.diff(simplified, axis=0)
    theta = geo_angle(simplified)
    # Select the index of the points with the greatest theta
    # Large theta is associated with greatest change in direction.
    idx = np.where(theta>np.rad2deg(min_angle))[0]+1
    print("number of significant turns = %d"%len(idx))

    return points, simplified, idx

def cart_find_corners(lat_lon_file,tolerance=1,min_angle=np.pi*0.22):
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

    # Use the Ramer-Douglas-Peucker algorithm to simplify the path
    # http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm
    # Python implementation: https://github.com/sebleier/RDP/
#    import pdb; pdb.set_trace()
    cart_points,_ = latlon2cart(points[:,0],points[:,1],np.zeros([3,3]))
    simplified = np.array(rdp.rdp(np.array(cart_points).T.tolist(), tolerance))

    print("simplified number of points = %d"%len(simplified))

    # compute the direction vectors on the simplified curve
    directions = np.diff(simplified, axis=0)
    theta = cart_angle(directions)
    # Select the index of the points with the greatest theta
    # Large theta is associated with greatest change in direction.
    idx = np.where(theta>min_angle)[0]+1
    print("number of significant turns = %d"%len(idx))

    simplified,_ = cart2latlon(simplified[:,0],simplified[:,1],simplified[:,2],np.zeros([3,3]))

    return points, np.array(simplified).T, idx

def create_geo_gc_dist(geoid):
    """
    Creates a function which calculates the distance between a great circle defined by two points and a third point
    based on the provided geoid which must have methods Inverse and ArcDirect

    Parameters
    ----------
    geoid : Geodesic
        geoid surface on which to calculate distance must be some form of general ellipsoid or sphere

    Returns
    -------
    geo_gc_dist : func
        function which calculates distance between gc and third point
    """

    def geo_gc_dist(point, start, end, verbose=False):

        if np.all(np.equal(start, end)): #null case where start and end are a single point
            return geoid.Inverse(*point,*start)["a12"]

        gc = geoid.Inverse(*start,*end) #find gc to calc distance from
        s_gcp = geoid.ArcDirect(*start,gc["azi1"]-90,90) #find pole of that gc
        e_gcp = geoid.ArcDirect(*end,gc["azi2"]-90,90) #find pole of that gc

        if np.sign(s_gcp["lat2"])!=np.sign(e_gcp["lat2"]): #This never should happen but in case...
            e_gcp = geoid.ArcDirect(*end,gc["azi2"]+90,90) #find pole of that gc

        p_gcp = geoid.Inverse(*point,s_gcp["lat2"],s_gcp["lon2"])
        if p_gcp["azi2"]>max([s_gcp["azi2"],e_gcp["azi2"]]) or p_gcp["azi2"]<min([s_gcp["azi2"],e_gcp["azi2"]]): #case of outside endpoints
            return min([geoid.Inverse(*start,*point)["a12"],geoid.Inverse(*point,*end)["a12"]])
        else: #case of within end points
            return np.abs(90-p_gcp["a12"]) #find distance between pole and point and difference from pole to gc

    return geo_gc_dist

def cart_angle(dir):
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

def geo_angle(points,geoid=Geodesic(6371.,0.)):
    """
    Returns the angles between two great circles for a set of points

    Parameters
    ----------
    points : 2D-Ndarray
        lat-lon array of points on a sphere between which to calculate angles

    Returns
    ----------
    angles : 1D-Ndarray
        angles between the great circles defined by points of shape (N-1,), with each value between 0 and 180 degrees.
    """
    return np.abs(np.diff(list(map(lambda a,b,c,d: geoid.Inverse(a,b,c,d)["azi1"],points[:-1,0],points[:-1,1],points[1:,0],points[1:,1]))))

def kmg9_rdp(A, tol, dist, start_index=0, last_index=None, lsz=None):
    """
    Temperary function to create animation of RDP if necessary for heuristic purposes.
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import pyskew.plot_gravity as pgrav
    from time import time
    from rasterio.enums import Resampling
    import matplotlib.patches as patches
    down_sample_factor = 10
    window = [247.5,262.5,10.,25.]
    sandwell_files_path = "/home/kevin/Projects/PySkew/raw_data/gravity/Sandwell/"
    verbose = False
    x,y = np.array(A).T

    stk = []
    if isinstance(last_index,type(None)): last_index = len(A)-1
    stk.append([start_index, last_index])
    global_start_index = start_index
    indices = np.ones(last_index - start_index + 1, dtype=bool)
    it = 0

    while stk:

        start_index, last_index = stk.pop()

        fig = plt.figure()
        llcrnrlon=min(utl.convert_to_0_360(y[start_index:last_index]))-2 if min(utl.convert_to_0_360(y[start_index:last_index]))-2>=0 else 0
        llcrnrlat=min(x[start_index:last_index])-2 if min(x[start_index:last_index])-2>=-89 else -89
        urcrnrlon=max(utl.convert_to_0_360(y[start_index:last_index]))+2 if max(utl.convert_to_0_360(y[start_index:last_index]))+2<=359 else 359
        urcrnrlat=max(x[start_index:last_index])+2 if max(x[start_index:last_index])+2<=90 else 89
        ax = pg.create_basic_map(projection='merc', resolution="10m",llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon, center_lon=180., fig=fig)
        ax.plot(y, x, 'k-', label='original path',zorder=1,transform=ccrs.PlateCarree())
        sidx = np.argsort(lsz[:,1])
        ax.plot(lsz[sidx,0], lsz[sidx,1], color="tab:grey", linestyle='-', label='EPR',zorder=1001,transform=ccrs.PlateCarree(), linewidth=3)

        dmax = 0.0
        index = start_index

        for i in range(index + 1, last_index):
            if indices[i - global_start_index]:
                d = dist(A[i], A[start_index], A[last_index], verbose=verbose)
                if d > dmax:
                    index = i
                    dmax = d

        ax.plot([y[start_index],y[last_index]],[x[start_index],x[last_index]], color="cyan", linestyle='-', label='Test Segment', marker='o', mec="k", zorder=3,transform=ccrs.PlateCarree())
        ax.plot(y[indices], x[indices], color="tab:orange", linestyle=':', label='simplified path',zorder=2,transform=ccrs.PlateCarree())
        ax.scatter(y[index], x[index], color="r", marker="o", s = 30, label='turning points',zorder=4,transform=ccrs.PlateCarree(),edgecolor="k")
        print(index,dmax,A[index])

        print("Grav Window and Down sample: ",window,down_sample_factor)
        all_lons,all_lats,all_grav = pgrav.get_sandwell(window,down_sample_factor,resample_method=Resampling.average,sandwell_files_path=os.path.join(sandwell_files_path,"*.tiff"))

        print("Plotting Gravity")
        start_time = time()
        print("Grid Sizes: ",all_lons.shape,all_lats.shape,all_grav.shape)
        fcm = ax.contourf(all_lons, all_lats, all_grav, 60, cmap="Blues_r", alpha=.75, transform=ccrs.PlateCarree(), zorder=0, vmin=0, vmax=255)
        print("Runtime: ",time()-start_time)

        sidx = np.argsort(lsz[:,1])
        ax.plot(lsz[sidx,0], lsz[sidx,1], color="tab:grey", linestyle='-', label='EPR',zorder=1001,transform=ccrs.PlateCarree(), linewidth=3)

        ax.set_extent(window, ccrs.PlateCarree())

        if dmax > tol:
            stk.append([start_index, index])
            stk.append([index, last_index])
        else:
            for i in range(start_index + 1, last_index):
                indices[i - global_start_index] = False

        fig.savefig("/home/kevin/Projects/RDP/img/animation/%s.png"%("0"*(3-len(str(it)))+str(it)))
        plt.show()
        it += 1

    return indices


