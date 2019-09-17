import os,glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
#from mpl_toolkits.basemap import Basemap
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
from .utilities import *

def plot_chron_info(chrons_info, m, coord_0_360=False, chron_dir=os.path.join("..","raw_data","chrons","cande"), barckhausen_path=os.path.join("..",'raw_data','chrons','Barckhausen2013','GSFML.Barckhausen++_2013_MGR.picks.gmt'),**kwargs):
    #Create Chron markers
    for chron_info in chrons_info:
        chron,chron_color = chron_info
        chron_path = os.path.join(chron_dir,"cande.%s"%str(chron))
        if not os.path.isfile(chron_path): print("no file %s, so there is probably no chron %s in cande's model."%(chron_path,str(chron))); continue
        fchron = open(chron_path,'r')
        lines = fchron.readlines()
        fchron.close()
        entries=[[],[]]
        for line in lines[1:]:
            entry = line.split()
            if entry[0]=='>':
                if len(entries[0])<2 or len(entries[1])<2: entries=[[],[]]; continue
                lats = entries[1]
                lons = entries[0]
#                XYM = [None,None]
#                XYM[0],XYM[1] = m(lons,lats)
#                XYM[0]=remove_off_map_points(XYM[0])
#                XYM[1]=remove_off_map_points(XYM[1])
#                if len(XYM[0])!=len(XYM[1]): print("error plotting chron %s on map, skipping this sz"%str(chron)); continue
                m.plot(lons,lats,color=chron_color,transform=ccrs.Geodetic(),**kwargs)
                entries=[[],[]]
            else:
                if coord_0_360:
                    entries[0].append(float(entry[0])); entries[1].append(float(entry[1]))
                else:
                    entries[0].append(convert_to_180_180(entry[0])); entries[1].append(float(entry[1]))
        ccz,gcz = get_barckhausen_2013_chrons(barckhausen_path=barckhausen_path)
        if str(chron) in ccz.keys():
            if ccz[str(chron)]:
                lonlats = np.array(ccz[str(chron)])
                if coord_0_360: lonlats[:,0] = convert_to_0_360(lonlats[:,0])
#                XYM = [None,None]
#                XYM[0],XYM[1] = m(lonlats[:,0],lonlats[:,1])
#                XYM[0]=remove_off_map_points(XYM[0])
#                XYM[1]=remove_off_map_points(XYM[1])
#                if len(XYM[0])!=len(XYM[1]): print("error plotting chron %s on map, skipping this sz"%str(chron)); continue
                m.plot(lons,lats,color=chron_color,transform=ccrs.Geodetic(),**kwargs)
        if str(chron) in gcz.keys():
            if gcz[str(chron)]:
                lonlats = np.array(gcz[str(chron)])
                if coord_0_360: lonlats[:,0] = convert_to_0_360(lonlats[:,0])
#                XYM = [None,None]
#                XYM[0],XYM[1] = m(lonlats[:,0],lonlats[:,1])
#                XYM[0]=remove_off_map_points(XYM[0])
#                XYM[1]=remove_off_map_points(XYM[1])
#                if len(XYM[0])!=len(XYM[1]): print("error plotting chron %s on map, skipping this sz"%str(chron)); continue
                m.plot(lons,lats,color=chron_color,transform=ccrs.Geodetic(),**kwargs)

def create_basic_map(projection='ortho',resolution='110m',center_lat=0, center_lon=0, min_lat=-80, max_lat=80, stereo_bound_lat=60, label_grid=True, grid_spacing=10, fig=None, ax_pos=111, landcolor="grey", return_all=False,llcrnrlat=-80,urcrnrlat=80,llcrnrlon=-180,urcrnrlon=180):
    #Create Map
    if fig==None: fig = plt.gcf()
    if projection=='ortho':
        proj = ccrs.Orthographic(central_latitude=center_lat,central_longitude=center_lon)
        m = fig.add_subplot(ax_pos,projection=proj)
    elif projection=='merc':
        proj = ccrs.Mercator(central_longitude=center_lon,min_latitude=min_lat, max_latitude=max_lat)
        m = fig.add_subplot(ax_pos,projection=proj)
        m.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], ccrs.PlateCarree())
    elif projection=="moll":
        proj = ccrs.Mollweide(central_longitude=center_lon)
        m = fig.add_subplot(ax_pos,projection=proj)
    elif projection=='npstere':
        proj = ccrs.NorthPolarStereo(central_longitude=center_lon)
        m = fig.add_subplot(ax_pos,projection=proj)
        m.set_extent([-180, 180, 90, stereo_bound_lat], ccrs.PlateCarree())
        m = make_circular_ax(m)
    elif projection=='spstere':
        proj = ccrs.SouthPolarStereo(central_longitude=center_lon)
        m = fig.add_subplot(ax_pos,projection=proj)
        m.set_extent([-180, 180, -90, -stereo_bound_lat], ccrs.PlateCarree())
        m = make_circular_ax(m)
    elif projection=='stere':
        proj = ccrs.Stereographic(central_longitude=center_lon,central_latitude=center_lat)
        m = fig.add_subplot(ax_pos,projection=proj)
        m.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], ccrs.PlateCarree())
    else: print("the create_basic_map function is very basic and may need to be expanded to handle this projection"); return
    land = cfeature.NaturalEarthFeature('physical', 'land', resolution, edgecolor="black", facecolor=landcolor)
    m.add_feature(land)
    if label_grid and projection!="merc": label_grid=False
    gl = m.gridlines(crs=ccrs.PlateCarree(), draw_labels=label_grid, color="k", linestyle=":")
    gl.xlabels_top = False
    if return_all: return m,gl,proj,fig
    else: return m

def make_circular_ax(m):
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    m.set_boundary(circle, transform=m.transAxes)

    return m

def add_chron_info_to_legend(chrons_info, handles):
    for chron_info in chrons_info:
        chron,chron_color = chron_info
        handle = mlines.Line2D([], [], color=chron_color, label='chron %s'%str(chron))
        handles.append(handle)
    return handles

def plot_transit_locations(chron_to_analyse, chrons_info, results_directory, data_directory='.', lon_0=180, lat_0=0, annotate=True):

    chron,chron_color = chron_to_analyse
    chron_name = "chron%s"%(str(chron))

    cut_tracks_path=os.path.join(data_directory,"usable_tracks_and_intersects_for_%s.txt"%str(chron_name))
    if os.path.isfile(cut_tracks_path):
        transiting_locations = [(transiting_cut.split('\t')[2],'aero') if 'aero' in transiting_cut.split('\t')[0] else (transiting_cut.split('\t')[2],'ship') for transiting_cut in open(cut_tracks_path,'r').readlines()]

    #Create Figure
    fig = plt.figure(figsize=(9, 9), dpi=80)

    #Create Map
    crossed_anom_map = create_basic_map(center_lat=lat_0,center_lon=lon_0)
#    crossed_anom_map = create_basic_map(projection='merc',resolution='l',llcrnrlon=360-120, llcrnrlat=-35, urcrnrlon=360-70, urcrnrlat=30)
#    crossed_anom_map.drawparallels([0],labels=[1,0,0,0])
#    crossed_anom_map.drawparallels([-10,10],color='b',labels=[1,0,0,0])
#    crossed_anom_map.drawparallels([-20,20],color='r',labels=[1,0,0,0])
#    crossed_anom_map.drawmeridians(np.arange(360-130,360-60,10),linewidth=0.1,labels=[0,0,0,1])

    #Create Chron markers
    plot_chron_info(chrons_info,crossed_anom_map,coord_0_360=True)

    if annotate:
        plt.annotate('Surveyor', xy=crossed_anom_map(360-150,46))
        plt.annotate('Mendocino', xy=crossed_anom_map(360-180,37))
        plt.annotate('Murray', xy=crossed_anom_map(360-160,29))
        plt.annotate('Molokai', xy=crossed_anom_map(360-155,19))
        plt.annotate('Clarion', xy=crossed_anom_map(360-140,14))
        plt.annotate('Cliperton', xy=crossed_anom_map(360-140,5))
        plt.annotate('Galapagos', xy=crossed_anom_map(360-150,-3))
        plt.annotate('Marquesas', xy=crossed_anom_map(360-150,-14))
        plt.annotate('Easter', xy=crossed_anom_map(360-150,-25))

    XYM = [None,None]
    handles,handle_types_saved = [],[]
    for loc in transiting_locations:
        if loc[1]=='aero': marker,edgecolor,facecolor,label='x','black','black',"Aeromag Data Crossing"
        else:  marker,edgecolor,facecolor,label='+','grey','grey',"Shipmag Data Crossing"
        anom_inter_lon = read_idx_from_string(loc[0])[0][1][0]
        anom_inter_lat = read_idx_from_string(loc[0])[0][1][1]
#        XYM[0],XYM[1] = crossed_anom_map(anom_inter_lon,anom_inter_lat)
        data_handle = crossed_anom_map.scatter(anom_inter_lon,anom_inter_lat,facecolors=facecolor,edgecolors=edgecolor,marker=marker,zorder=3,label=label,transform=ccrs.Geodetic())
        if loc[1]=='aero' and 'aero' not in handle_types_saved:
            handles.append(data_handle); handle_types_saved.append('aero')
        elif loc[1]=='ship' and 'ship' not in handle_types_saved:
            handles.append(data_handle); handle_types_saved.append('ship')

    plt.title("Transited Spreading Zones")
    add_chron_info_to_legend(chrons_info,handles)
    legend = plt.legend(handles=handles,loc='best')
    frame = legend.get_frame()
    frame.set_alpha(.6)
    fig.savefig(os.path.join(results_directory,"transited_spreading_zones.png"))
    plt.close(fig)

def plot_tracks(chrons_info, results_directory, tracks=[], track_dir="all_tracks",lon_0=180,lat_0=0,cuts=False):
    """
    plots track files in tracks with chron info on a default ortho map centered on the pacific for rough estimation of intersept and places these plots in the results directory. If no tracks provided it plots them all
    """

    if tracks==[]: tracks = glob.glob('../raw_data/hi_alt/**/*.DAT')+glob.glob('../raw_data/ship/**/*.lp')

    all_tracks_dir=os.path.join(results_directory,track_dir)
    check_dir(all_tracks_dir)

    #Start loop making plots
    for track in tracks:

        #Get Track Name
        track_name = os.path.basename(os.path.basename(track))
        if not cuts: track_name = track_name.split('.')[0]
        else: track_name = track_name.replace('.','-')

        #Create Figure
        fig = plt.figure(figsize=(9, 9), dpi=80)

        #Create Map
        aero_track_map = create_basic_map(projection='ortho',center_lon=lon_0,center_lat=lat_0,fig=fig)

        #Create Chron markers
        plot_chron_info(chrons_info,aero_track_map)

        dfin = open_mag_file(track)

        lats = list(map(float,dfin['lat'].tolist()))
        lons = list(map(float,dfin['lon'].tolist()))

#        XYM = [None,None]
#        XYM[0],XYM[1] = aero_track_map(lons,lats)
#        XYM[0] = remove_off_map_points(XYM[0])
#        XYM[1] = remove_off_map_points(XYM[1])
#        if len(XYM[0])==0 or len(XYM[1])==0: continue
        aero_track_handle, = aero_track_map.plot(lons,lats,color='k',zorder=3,label=track_name,transform=ccrs.Geodetic())

        #plot title and labels
        plt.title(track_name)
        handles=[aero_track_handle]
        add_chron_info_to_legend(chrons_info, handles)
        plt.legend(handles=handles,loc='best')

        fig.savefig(os.path.join(all_tracks_dir,track_name+".png"))
        plt.close(fig)

def plot_track_cuts(x,y,sx,sy,idx,chrons_info,track_name,results_directory):
    fig = plt.figure(figsize=(9, 9), dpi=80) #(comment from here down to stop plotting, so it runs faster)

    llcrnrlon=min(convert_to_0_360(y))-20 if min(convert_to_0_360(y))-20>=0 else 0
    llcrnrlat=min(x)-20 if min(x)-20>=-89 else -89
    urcrnrlon=max(convert_to_0_360(y))+20 if max(convert_to_0_360(y))+20<=359 else 359
    urcrnrlat=max(x)+20 if max(x)+20<=90 else 89

    #Make a map
    turning_points_map = create_basic_map(projection='merc',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon, fig=fig)

#            turning_points_map = create_basic_map()

    #Create Chron markers
    plot_chron_info(chrons_info,turning_points_map,coord_0_360=True)

    #translate lat lon data into meters (I think, could just be random map values)
#    XYM = [None,None]
#    Y,X = turning_points_map(convert_to_0_360(y),x)
#    SY,SX = turning_points_map(convert_to_0_360(sy),sx)
#    X=np.array(remove_off_map_points(X))
#    Y=np.array(remove_off_map_points(Y))
#    SX=np.array(remove_off_map_points(SX))
#    SY=np.array(remove_off_map_points(SY))

    #Plot the data
    if len(x)!=len(y) or len(sx)!=len(sy): print("error plotting track on map it probably crosses a pole, opening debugger"); import pdb; pdb.set_trace()
    handle1, = turning_points_map.plot(y, x, 'b-', label='original path',zorder=1,transform=ccrs.Geodetic())
    handle2, = turning_points_map.plot(sy, sx, 'g--', label='simplified path',zorder=2,transform=ccrs.Geodetic())
    handle3, = turning_points_map.plot(sy[idx], sx[idx], 'ro', markersize = 7, label='turning points',zorder=3,transform=ccrs.Geodetic())

    #Name stuff
    plt.title(track_name)
    handles=[handle1,handle2,handle3]
    add_chron_info_to_legend(chrons_info,handles)
    plt.legend(handles=handles,loc='best')

    #Save the plot
    fig.savefig(os.path.join(results_directory,"turning_points",track_name+".png"))
    plt.close(fig)

def plot_az_strike(track,spreading_zone_file,idx,az,strike,chron_color,chron_name,results_directory,fout_name):

    #Create Figure
    fig = plt.figure(figsize=(9, 9), dpi=80)

    #Create Chron markers
    dft = open_mag_file(track)
    lt = [[convert_to_0_360(lon),float(lat)] for lon,lat in zip(dft['lon'],dft['lat'])]
    at = np.array(lt)
    lsz = [list(map(float,line.split())) for line in open(spreading_zone_file).readlines()]
    asz = np.array(lsz)

    #Create Map
#            gcm = create_basic_map() #uses defaults, hit shift-tab in parens to see what they are

    llcrnrlon=min(at[:,0])-20 if min(at[:,0])-20>0 else 0
    llcrnrlat=min(at[:,1])-20 if min(at[:,1])-20>-89 else -89
    urcrnrlon=max(at[:,0])+20 if max(at[:,0])+20<360 else 360
    urcrnrlat=max(at[:,1])+20 if max(at[:,1])+20<89 else 89

    gcm = create_basic_map(projection='merc',llcrnrlat=llcrnrlat,urcrnrlat=urcrnrlat,llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon, fig=fig)

    sz_handle, = gcm.plot(asz[:,0],asz[:,1],color=chron_color,zorder=1,label=chron_name,transform=ccrs.Geodetic())

    gcm_handle, = gcm.plot(at[:,0],at[:,1],color='k',zorder=2,label=os.path.basename(track),transform=ccrs.Geodetic())

    gcm.scatter(at[idx[1][0]][0],at[idx[1][0]][1],color='g',marker='o',s=10,zorder=3,label='nearest intercept',transform=ccrs.Geodetic())

    geodict = Geodesic.WGS84.Direct(float(at[idx[1][0]][1]), float(at[idx[1][0]][0]), float(az), 1000000)
    b_lon,b_lat = (360+geodict["lon2"])%360, geodict["lat2"]
    gcm.arrow(b_lon,b_lat, b_lon-at[idx[1][0]][0], b_lat-at[idx[1][0]][1], fc="white", ec="r", linewidth=1, head_width=1, head_length=1, label='azimuth',transform=ccrs.Geodetic())

    geodict = Geodesic.WGS84.Direct(float(at[idx[1][0]][1]), float(at[idx[1][0]][0]), float(strike), 1000000)
    b_lon,b_lat = (360+geodict["lon2"])%360, geodict["lat2"]
    gcm.arrow(b_lon,b_lat, b_lon-at[idx[1][0]][0], b_lat-at[idx[1][0]][1], fc="white", ec="pink", linewidth=1, head_width=1, head_length=1, label='strike',transform=ccrs.Geodetic())

    #plot title and labels
    plt.title(os.path.basename(track))
    plt.legend(loc='best')

    az_plots_dir = os.path.join(results_directory,"azimuth_strike_plots")
    check_dir(az_plots_dir)

    fig.savefig(os.path.join(az_plots_dir,os.path.basename(fout_name)[:-5]+"png"))
    plt.close(fig)
