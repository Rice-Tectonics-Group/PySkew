import os,sys
import subprocess
import pandas as pd
import numpy as np
import pmagpy.ipmag as ipmag
import matplotlib.pyplot as plt
if 'darwin' in sys.platform:
    try:
        plt.switch_backend('QT5Agg')
    except ImportError as e:
        print("QT5 not found, this library is optional on MacOS as it allows for better data visualization moving on with default graphics system.")
import matplotlib.lines as mlines
from collections import OrderedDict
from matplotlib.patches import Rectangle,ConnectionPatch,Polygon
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
from .skewness import *
from .plot_geographic import *
from .utilities import *

def plot_all_lunes_seperate(deskew_path,transform=ccrs.Geodetic()):

    comps = filter_deskew_and_calc_aei(deskew_path)
    results_dir = comps["results_dir"].iloc[0]

    lunes_dir = os.path.join(results_dir,"lunes_plots")
    check_dir(lunes_dir)

    if "E" in comps["comp_name"].iloc[0]:
        first_comp_dir = "E"
        last_comp_dir = "V"
    elif "V" in comps["comp_name"].iloc[0]:
        first_comp_dir = "V"
        last_comp_dir = "E"
    else: print("invalid name for the first component %s, check data. Aborting"%comps["comp_name"]); return

    # For every crossing...
    for i,row in comps.iterrows():
        track_name = row["comp_name"].rstrip(".Ed .Vd .lp").replace(".","-")

        if row["track_type"]=='ship' or first_comp_dir in row["comp_name"]:
            # Create figure
            fig = plt.figure(figsize=(16,9), dpi=200)

            # Create map
            gcm,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0,fig=fig,return_all=True)

        # Find array of great semicircle azimuths (degrees)
        azi = np.linspace(float(row["strike"])-(180),float(row["strike"]),100)

        # For first bounding azimuth...
        # Find inclination (degrees)
        inc = np.degrees(np.arctan(np.tan(np.deg2rad(row["aei"]))*np.sin(np.deg2rad(azi+180-float(row["strike"])))))
        # Find paleocolatitude (degrees)
        clt = 90-np.degrees(np.arctan(np.tan(np.deg2rad(inc))*0.5))
        # Find great circle endpoint
        gc_points_and_azis = [Geodesic.WGS84.ArcDirect(float(row["inter_lat"]),float(row["inter_lon"]), azi1, clt1) for clt1,azi1 in zip(clt,azi)]
        gc_lon = [gcd["lon2"] for gcd in gc_points_and_azis]
        gc_lat = [gcd["lat2"] for gcd in gc_points_and_azis]

        if row["track_type"]=='aero':
            if "V" in row["comp_name"]: linestyle,label="--","Aero Vertical"
            else: linestyle,label="-","Aero East"
        elif row["track_type"]=='ship':
            linestyle,label=":","Ship"

        # Draw great circle
        gcm.plot(gc_lon, gc_lat, color='black',linestyle=linestyle,label=label,transform=transform)

        if row["track_type"]=='ship' or last_comp_dir in row["comp_name"]:
            plt.title(track_name)
            plt.legend(loc='best')

            plot_path = os.path.join(lunes_dir,track_name+".png")

            print("saving to %s"%plot_path)
            fig.savefig(plot_path)
            plt.close(fig)

def plot_gc_from_pole(slat,slon,plat,plon,length=360,error=None,m=None,spacing=1,**kwargs):
    geodict = Geodesic.WGS84.Inverse(slat,slon,plat,plon)
    azi_mean = geodict["azi1"]
    if isinstance(error,int) or isinstance(error,float):
        geodict = Geodesic.WGS84.ArcDirect(plat,plon,geodict['azi2']-90,error)
        geodict = Geodesic.WGS84.Inverse(slat,slon,geodict['lat2'],geodict['lon2'])
        azi_a95 = (azi_mean-geodict['azi1'])/np.sqrt(2)
    else: azi_a95 = None

    m = plot_great_circle(slon,slat,azi_mean,error=azi_a95,length=length,m=m,spacing=spacing,**kwargs)
    return m

def plot_great_circle(plon,plat,azimuth,length=360,error=None,m=None,geoid=Geodesic.WGS84,spacing=1,transform=ccrs.Geodetic(),**kwargs):

    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=200)

        # Create map
        m,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0, fig=fig, return_all=True)

    lats,lons = [],[]
    for dis in np.arange(0,length,spacing):
        geodict = geoid.ArcDirect(plat,plon,azimuth,dis)
        lons.append(geodict["lon2"])
        lats.append(geodict["lat2"])
    m.plot(lons,lats,transform=transform,**kwargs)

    if isinstance(error,int) or isinstance(error,float):
        if "linestyle" in kwargs.keys(): kwargs.pop("linestyle")
        lats,lons = [],[]
        for dis in np.arange(0,length,spacing):
            geodict = geoid.ArcDirect(plat,plon,azimuth+error,dis)
            lons.append(geodict["lon2"])
            lats.append(geodict["lat2"])
        m.plot(lons,lats,transform=transform,linestyle='--',**kwargs)

        lats,lons = [],[]
        for dis in np.arange(0,length,spacing):
            geodict = geoid.ArcDirect(plat,plon,azimuth-error,dis)
            lons.append(geodict["lon2"])
            lats.append(geodict["lat2"])
        m.plot(lons,lats,transform=transform,linestyle='--',**kwargs)

    return m

def plot_small_circles(plon,plat,m=None,range_arcdis=(10,180,10),range_azis=(-180,180,1),**kwargs):
    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=200)
        # Create map
        m,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0, fig=fig, return_all=True)

    for arcdis in range(*range_arcdis):
        plot_small_circle(plon,plat,arcdis,m=m,range_azis=range_azis,**kwargs)

    return m

def plot_small_circle(plon,plat,arcdis,error=None,m=None,range_azis=(-180,180,1), alpha_line=1.0, filled=True, geoid=Geodesic(6731.,0.), transform=ccrs.Geodetic(), **kwargs):
    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=200)
        # Create map
        m,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0, fig=fig, return_all=True)

    sml_circ_points = []
    for azi in np.arange(*range_azis):
        geo_dict = geoid.ArcDirect(plat,plon,azi,arcdis)
        sml_circ_points.append([utl.convert_to_0_360(geo_dict['lon2']),geo_dict['lat2']])
    sml_circ_points.append(sml_circ_points[0])

    if isinstance(error,float) or isinstance(error,int):
        lb_sml_circ_points,ub_sml_circ_points = [],[]
        for azi in range(*range_azis):
            geo_dict = geoid.ArcDirect(plat,plon,azi,arcdis+error)
            ub_sml_circ_points.append([utl.convert_to_0_360(geo_dict['lon2']),geo_dict['lat2']])
            geo_dict = geoid.ArcDirect(plat,plon,azi,arcdis-error)
            lb_sml_circ_points.append([utl.convert_to_0_360(geo_dict['lon2']),geo_dict['lat2']])
        ub_sml_circ_points.append(ub_sml_circ_points[0])
        lb_sml_circ_points.append(lb_sml_circ_points[0])

        ub_asml_circ_points = np.array(ub_sml_circ_points)

        lb_asml_circ_points = np.array(lb_sml_circ_points)

        if filled:
            if 'linewidth' in kwargs.keys():
                linewidth = kwargs['linewidth']
                kwargs.pop('linewidth')
            m.fill_between(ub_asml_circ_points[:,0],lb_asml_circ_points[:,1],ub_asml_circ_points[:,1],transform=transform,linewidth=0,**kwargs)
            kwargs['linewidth'] = linewidth
        else:
            m.plot(ub_asml_circ_points[:,0],ub_asml_circ_points[:,1],transform=transform,linestyle='--',**kwargs)
            m.plot(lb_asml_circ_points[:,0],lb_asml_circ_points[:,1],transform=transform,linestyle='--',**kwargs)

    asml_circ_points = np.array(sml_circ_points)
    if "alpha" in kwargs.keys(): kwargs.pop("alpha")
    m.plot(asml_circ_points[:,0][:-1],asml_circ_points[:,1][:-1],transform=transform, alpha=alpha_line, **kwargs)

    return m

def plot_pole_track_from_files(ellipse_paths,m=None,**kwargs):
    ellipse_data = []
    for ellipse_path in ellipse_paths:
        lon,lat,az,a,b = open_ellipse_file(ellipse_path)
        ellipse_data.append([lon,lat,az,a,b])
    m,_ = plot_pole_track(ellipse_data,m=m,**kwargs)
    return m

def get_alpha_func(max_error,min_error,max_alpha,min_alpha):
    A = (max_alpha-min_alpha)/(min_error-max_error)
    B = max_alpha-A*min_error
    calc_alpha = lambda x : A*x+B if (A*x+B)>0 else min_alpha
    return calc_alpha

def get_alpha_func_from_ellipse_data(ellipse_data,max_alpha,min_alpha):
    error_array = abs(np.array(ellipse_data)[:,3])+abs(np.array(ellipse_data)[:,4])
    max_error,min_error = max(error_array),min(error_array)
    calc_alpha = get_alpha_func(max_error,min_error,max_alpha,min_alpha)
    return calc_alpha

def plot_north_pole(m):
    m.scatter(0, 90, facecolor='black', edgecolor='black', marker='+', s=30, zorder=120)
    return m

def plot_pole(lon,lat,az,a,b,m=None,color='cyan',zorder=3,pole_text=None,pole_text_pos=None,alpha_all=False,filled=True,transform=ccrs.PlateCarree(),**kwargs):
    #a and b should be full axes of the ellipse not semi-axes

    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=200)
        # Create map
        m,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0, fig=fig, return_all=True)

    if "edgecolors" not in kwargs.keys(): kwargs["edgecolors"]="black"
    if "facecolors" not in kwargs.keys(): kwargs["facecolors"]=color
    if kwargs.get("edgecolors") == "black" and kwargs.get("facecolors") != color: color = kwargs.get("facecolors")
    elif kwargs.get("facecolors") == "white" and kwargs.get("edgecolors") != color: color = kwargs.get("edgecolors")
    if "alpha" in kwargs.keys(): alpha = kwargs.pop("alpha")
    else: alpha = None

    if alpha_all: m.scatter(lon, lat,  zorder=zorder,transform=transform, alpha=alpha, **kwargs)
    else: m.scatter(lon, lat,  zorder=zorder,transform=transform, **kwargs)
    if pole_text!=None:
        if pole_text_pos!=None: plt.text(*m(*pole_text_pos),pole_text,zorder=500,transform=transform)
        else:
            tazi = 90
            geodict=Geodesic.WGS84.ArcDirect(lat,lon,tazi,1)
            plt.text(geodict["lon2"],geodict["lat2"],pole_text,zorder=500,va='top',ha='center',transform=transform)
    if "edgecolors" in kwargs: kwargs.pop("edgecolors")
    if "facecolors" in kwargs: kwargs.pop("facecolors")
    if "marker" in kwargs: kwargs.pop("marker")
    if "linewidths" in kwargs: kwargs.pop("linewidths")
    if "s" in kwargs: kwargs.pop("s")
    if filled: ipmag.ellipse(m, lon, lat, (a*111.11), (b*111.11), az, n=360, filled=filled, facecolor=color, edgecolor='#00000088', zorder=zorder-1, alpha=alpha, **kwargs)
    else: ipmag.ellipse(m, lon, lat, (a*111.11), (b*111.11), az, n=360, filled=filled, color=color, zorder=zorder-1, alpha=alpha, **kwargs)

    return m

def plot_apw_legend(m,**kwargs):
    handles,labels = m.get_legend_handles_labels()
    handles += [mlines.Line2D([0],[0],marker='o', markeredgecolor='k', markerfacecolor='grey', label="Normal Polarity",linestyle=''),
                            mlines.Line2D([0],[0],marker='s', markeredgecolor='k', markerfacecolor='grey', label='Mixed Polarity',linestyle=''),
                            mlines.Line2D([0],[0],marker='o', markeredgecolor='grey', markerfacecolor='w', label="Reversed Polarity",linestyle=''),
                            mlines.Line2D([0],[0],marker='^', markeredgecolor='k', markerfacecolor='grey', label='Mixed Polarity',linestyle=''),
                            mlines.Line2D([0],[0],marker='D', markeredgecolor='k', markerfacecolor='grey', label='Mixed Polarity',linestyle=''),
                            mlines.Line2D([0],[0],marker='v', markeredgecolor='k', markerfacecolor='grey', label='Mixed Polarity',linestyle='')]
    labels += ["Normal Polarity","Mixed Polarity","Reversed Polarity","Normal Colatitude/Azimuth","Mixed Colatitude/Azimuth","Reversed Colatitude/Azimuth"]
    by_label = OrderedDict(zip(labels, handles))
    leg = m.legend(by_label.values(), by_label.keys(),**kwargs)
    return m,leg

def plot_pole_track(ellipse_data, m=None,modify_alpha_with_error=True, min_alpha=0.2, max_alpha=0.8, annotations=[],annotation_positions=[], transform=ccrs.Geodetic(), **kwargs):

    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=200)
        # Create map
        m,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0, fig=fig, return_all=True)

    if modify_alpha_with_error:
        if len(ellipse_data)>1: calc_alpha = get_alpha_func_from_ellipse_data(ellipse_data,max_alpha,min_alpha)
        else: calc_alpha = lambda x: max_alpha

    lats,lons,label,zorder = [],[],None,4
    if 'label' in kwargs.keys(): label = kwargs.pop('label')
    if 'zorder' in kwargs: zorder = kwargs.pop('zorder')+2
    for i,(lon,lat,az,a,b) in enumerate(ellipse_data):
        if modify_alpha_with_error: kwargs['alpha'] = calc_alpha(a+b)
        if i < len(annotations): pole_text = annotations[i]
        else: pole_text = None
        if i < len(annotation_positions): pole_text_pos = annotation_positions[i]
        else: pole_text_pos = None
        m = plot_pole(lon,lat,az,a,b,m=m,zorder=zorder,pole_text=pole_text,pole_text_pos=pole_text_pos,**kwargs)
        lons.append(lon);lats.append(lat)

    m.plot(lons,lats,label=label,zorder=zorder-1,transform=transform,**kwargs)
    m.plot(lons,lats,color='k',zorder=zorder-2,linewidth=1.5,transform=transform)

    return m,plt.gcf()

def plot_pole_tracks(tracks, plot_kwargs=None, m=None,modify_alpha_with_error=True, min_alpha=0.2, max_alpha=0.8, return_alpha_calculator=False):

    if modify_alpha_with_error:
        error_array = abs(np.array([pole[3] for track in tracks for pole in track]))+abs(np.array([pole[3] for track in tracks for pole in track]))
        max_error,min_error = max(error_array),min(error_array)
        calc_alpha = get_alpha_func(max_error,min_error,max_alpha,min_alpha)
    else: calc_alpha = lambda x: max_alpha

    if plot_kwargs==None: plot_kwargs = [{} for i in range(len(tracks))]

    for track,kwargs in zip(tracks,plot_kwargs):
        if modify_alpha_with_error:
            error_array = abs(np.array(track)[:,3])+abs(np.array(track)[:,4])
            kwargs['min_alpha'],kwargs['max_alpha'] = calc_alpha(max(error_array)),calc_alpha(min(error_array))
        m,_ = plot_pole_track(track,m=m,modify_alpha_with_error=modify_alpha_with_error,**kwargs)

    if return_alpha_calculator and modify_alpha_with_error:
        return m,calc_alpha
    else:
        return m

def plot_pole_with_lunes(deskew_path,ellipse_path):

    comps = filter_deskew_and_calc_aei(deskew_path)

    # Create figure
    fig = plt.figure(figsize=(16,9), dpi=300)

    # Create map
    gcm,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0,fig=fig,return_all=True,resolution='10m')

    lon,lat,az,a,b = open_ellipse_file(ellipse_path)

    plot_lunes(comps,gcm,pole_lat=lat)

    plot_pole(lon,lat,az,a,b,m=gcm,zorder=5)

    out_path = os.path.join(comps['results_dir'].iloc[0],"poles.png")
    check_dir(os.path.dirname(out_path))
    print("saving to %s"%out_path)
    fig.savefig(out_path)
    plt.close(fig)

def plot_lunes_and_save(deskew_path):

    comps = filter_deskew_and_calc_aei(deskew_path)

    # Create figure
    fig = plt.figure(figsize=(16,9), dpi=200)

    # Create map
    gcm,gl,proj,fig = create_basic_map(projection='npstere', center_lon=0,fig=fig,return_all=True)

    plot_lunes(comps,gcm)

    out_path = os.path.join(comps['results_dir'].iloc[0],"lunes.png")
    print("saving to %s"%out_path)
    fig.savefig(out_path)
    plt.close(fig)


def plot_lunes(comps, gcm, idx_selected=None, average_vector_data=False, plot_legends=True, geoid=Geodesic(6371000.,0.), transform=ccrs.Geodetic(),**kwargs):

    # For every crossing...
    for i,row in comps.iterrows():
        if row["comp_name"].startswith('#'): continue

        aei,strike,lat,lon = float(row["aei"]),convert_to_0_360(row['strike']),float(row["inter_lat"]),float(row["inter_lon"])
        if row["track_type"]=='aero':
            if "Ed.lp" in row["comp_name"]:
                linestyle="-"
                if average_vector_data:
                    other_comp = row["comp_name"].replace(".Ed.lp",".Vd.lp")
                    oth_row = comps[comps["comp_name"]==other_comp].iloc[0]
                    lat = (row["inter_lat"] + oth_row["inter_lat"])/2
                    lon = (row["inter_lon"] + oth_row["inter_lon"])/2
                    strike = (row["strike"] + oth_row["strike"])/2
                    aei = (row["aei"] + oth_row["aei"])/2
            elif "Hd.lp" in row["comp_name"]:
                linestyle="-"
                if average_vector_data:
                    other_comp = row["comp_name"].replace(".Hd.lp",".Vd.lp")
                    oth_row = comps[comps["comp_name"]==other_comp].iloc[0]
                    lat = (row["inter_lat"] + oth_row["inter_lat"])/2
                    lon = (row["inter_lon"] + oth_row["inter_lon"])/2
                    strike = (row["strike"] + oth_row["strike"])/2
                    aei = (row["aei"] + oth_row["aei"])/2
            elif average_vector_data: continue
            else: linestyle="--"
            linewidth=1
        if row["track_type"]=='ship':
            linestyle,linewidth=":",1

        if aei<90 and aei>-90: azi = (360+np.linspace(strike-180,strike,180))%360 #Possible Declinations
        else: azi = (360+np.linspace(strike,strike+180,180))%360 #Possible Declinations

        # For first bounding azimuth...
        # Find inclination (degrees)
        inc = (np.arctan(np.tan(np.deg2rad(aei))*np.sin(np.deg2rad(azi+180-strike)))) #The 180 is because of the right hand rule and the sign to handle northern hemisphere
        # Find paleocolatitude (degrees)
        clt = 90-np.rad2deg(np.arctan2(np.tan(inc),2))
        # Find great circle points
        gc_points_and_azis = [geoid.ArcDirect(lat,lon, azi1, clt1) for clt1,azi1 in zip(clt,azi)]
        gc_lon = [gcd["lon2"] for gcd in gc_points_and_azis]
        gc_lat = [gcd["lat2"] for gcd in gc_points_and_azis]

        # Draw great semi-circle
        gcm.scatter([gc_lon[0],gc_lon[-1]], [gc_lat[0],gc_lat[-1]], edgecolor='k', facecolor='none', zorder=10,transform=ccrs.PlateCarree())
        if (not isinstance(idx_selected,type(None))) and i==idx_selected: color = "#FF6C6C"
        else: color = (float(row["r"]),float(row["g"]),float(row["b"]))
        gcm.plot(gc_lon, gc_lat, color=color, linestyle=linestyle, linewidth=linewidth,transform=transform,**kwargs)

    comps["inter_lat"] = comps["inter_lat"].apply(float)
    tmp_comps = comps.sort_values(by="inter_lat",ascending=False)

    if plot_legends:
        sz_handles=[]
        for i,row in tmp_comps[["sz_name","r","g","b"]].drop_duplicates().iterrows():
            sz_handle = mlines.Line2D([], [], color=(float(row['r']),float(row['g']),float(row['b'])), label=row['sz_name'])
            sz_handles.append(sz_handle)

        track_type_handles=[]
        for track_type,style in [['Aeromag East','-'],['Aeromag Vertical','--'],['Shipmag',':']]:
            track_type_handle = mlines.Line2D([], [], color='k', linestyle=style, label=track_type)
            track_type_handles.append(track_type_handle)

    #    handles, labels = plt.gca().get_legend_handles_labels()
    #    by_label = OrderedDict(zip(labels, handles))
    #    legend = plt.legend(by_label.values(), by_label.keys(),loc=1)

        track_type_legend = plt.legend(handles=track_type_handles,loc=2,framealpha=.7)
    #    sz_legend = plt.legend(handles=sz_handles,bbox_to_anchor=(.95, .91),bbox_transform=plt.gcf().transFigure, frameon=False)
        sz_legend = plt.legend(handles=sz_handles,loc=1,framealpha=.7)

        plt.gca().add_artist(track_type_legend)

    return gcm

def plot_chron_span_on_axes(sz_name, axes, anom_age_span,spreading_rate_path="../raw_data/spreading_rate_model.txt"):

    srf = generate_spreading_rate_model(spreading_rate_path)[0]
    step = .01
    srf_sz = lambda x: step*srf(sz_name,x)
    anom_width = sum(map(srf_sz,np.arange(*anom_age_span,step)))

    for axis in axes:
        axis.axvspan(-anom_width/2, anom_width/2, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)

def plot_synthetic(sz_name, anom_age_span, ax=plt.gca(),timescale_path="../raw_data/timescale_gradstein2012.txt",spreading_rate_path="../raw_data/spreading_rate_model.txt",age_min=40.0,age_max=100.0,layer_depth=4.5,layer_thickness=0.5,layer_mag=1000.,azi=90.,rd=0.,ri=90.,ad=0.,ai=90.,fix_sta=False,fix_end=False,twf=0.,length=4096, xlims=[-500,500], clip_on=False, **kwargs):
    mag_syn,dis_syn,samp_dis = make_synthetic(age_min,age_max,layer_depth,layer_thickness,layer_mag,azi,rd,ri,ad,ai,fix_sta,fix_end,twf,timescale_path,sz_name=sz_name,spreading_rate_path=spreading_rate_path,length=length)

    #Center Synthetic
    srf = generate_spreading_rate_model(spreading_rate_path)[0]
    if age_min!=age_max: step = (age_max-age_min)/(length-1)
    else: step = .01
    srf_sz = lambda x: step*srf(sz_name,x)
    center_dis = sum(map(srf_sz,np.arange(0,sum(anom_age_span)/2+step,step)))
#    dis_anom_min = sum(map(srf_sz,np.arange(0,age_min+step,step)))
#    dis_anom_max = sum(map(srf_sz,np.arange(0,age_max+step,step)))
#    neg_anom_max = sum(map(srf_sz,np.arange(-age_max+step,0,step)))
#    neg_anom_min = sum(map(srf_sz,np.arange(-age_min+step,0,step)))
#    anom_width = abs(dis_anom_max-dis_anom_min)/2
#    center_dis = dis_anom_max - anom_width
#    neg_anom_width = abs(neg_anom_max-neg_anom_min)/2
#    neg_anom = -(neg_anom_max - neg_anom_width) - center_dis
    dis_syn = np.array(dis_syn) - center_dis
    if not clip_on:
        left_idx = np.abs(dis_syn - xlims[0]).argmin()
        right_idx = np.abs(dis_syn - xlims[1]).argmin()
        dis_syn,mag_syn = dis_syn[left_idx:right_idx],mag_syn[left_idx:right_idx]

#    zln = ax.plot(-dis_syn,np.zeros(len(dis_syn)),'k--')
    sln = ax.plot(-dis_syn,mag_syn, clip_on=clip_on,**kwargs) #reverse because I'm dumb and changed convention

#    if not clip_on:
#        ax.set_xlim(xlims)
#        start_box = ax.transAxes.transform([0,0])[0]
#        end_box = ax.transAxes.transform([1,0])[0] - start_box
#        clip_patch = Rectangle([start_box,-1e7],end_box,1e9)
#        ax.lines[0].set_clip_path(clip_patch)
#        ax.lines[1].set_clip_path(clip_patch)

    return min(dis_syn),max(dis_syn)

def plot_scale_bars(ax,size_x=100,size_y=100,x_unit='km',y_unit='nT', scale_bar_figure_location = (.95,.05), offset_of_bars = .02):
    #do figure to data calculations to get the proper width and location of bars
    xbar_xy_figure_location = (scale_bar_figure_location[0]-offset_of_bars,scale_bar_figure_location[1])
    ybar_xy_figure_location = (scale_bar_figure_location[0],scale_bar_figure_location[1]+offset_of_bars)
    xbar_xy_display = ax.get_figure().transFigure.transform(xbar_xy_figure_location)
    ybar_xy_display = ax.get_figure().transFigure.transform(ybar_xy_figure_location)
    xbar_xy_data = ax.transData.inverted().transform(xbar_xy_display)
    ybar_xy_data = ax.transData.inverted().transform(ybar_xy_display)
    xbar_xytext_data = ((xbar_xy_data[0]-size_x),xbar_xy_data[1])
    xbar_centered_text_data = ((xbar_xytext_data[0]+size_x/2),xbar_xytext_data[1])
    ybar_xytext_data = (ybar_xy_data[0],(ybar_xy_data[1]+size_y))
    xbar_xytext_display = ax.transData.transform(xbar_xytext_data)
    xbar_centered_text_display = ax.transData.transform(xbar_centered_text_data)
    ybar_xytext_display = ax.transData.transform(ybar_xytext_data)
    xbar_xytext_figure_location = ax.get_figure().transFigure.inverted().transform(xbar_xytext_display)
    xbar_centered_text_figure_location = ax.get_figure().transFigure.inverted().transform(xbar_centered_text_display)
    ybar_xytext_figure_location = ax.get_figure().transFigure.inverted().transform(ybar_xytext_display)

    #create xbar
    ax.annotate('',xy=xbar_xy_figure_location, xycoords='figure fraction', xytext=(xbar_xytext_figure_location), textcoords='figure fraction', arrowprops=dict(arrowstyle='|-|'))
    ax.annotate('%d %s'%(size_x,x_unit),xy=xbar_xy_figure_location, xycoords='figure fraction', xytext=(xbar_centered_text_figure_location), textcoords='figure fraction', va='bottom', ha='center', fontsize=10)

    #create ybar
    ax.annotate('',xy=ybar_xy_figure_location, xycoords='figure fraction', xytext=(ybar_xytext_figure_location), textcoords='figure fraction', arrowprops=dict(arrowstyle='|-|'))
    ax.annotate('%d %s'%(size_y,y_unit),xy=ybar_xy_figure_location, xycoords='figure fraction', va='top', ha='center', fontsize=10)

def plot_skewness_data(deskew_row, phase_shift, ax, xlims=[-500,500], clip_on=False, return_objects=False, flip=False, **kwargs):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","dec_year","mag","lat","lon"],delim_whitespace=True)

    if flip:
        projected_distances = calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),(180+deskew_row['strike'])%360)
    else:
        projected_distances = calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),deskew_row['strike'])

    shifted_mag = phase_shift_data(data_df['mag'].tolist(),phase_shift)

    proj_dist = projected_distances['dist'].tolist()
    if deskew_row["quality"]=="g": #is "good" data
        if "color" in kwargs.keys(): pass
        else: kwargs["color"] = "k"
        zcolor = "k"
        zalpha = 1.0
    else: #is bad" data
        if "color" in kwargs.keys(): pass
        else: kwargs["color"] = "grey"
        zcolor = "grey"
        zalpha = .7
    zlinestyle = "--"
    if not clip_on and not any(list(map(lambda x: isinstance(x,type(None)),xlims))):
        left_idx = np.abs(np.array(proj_dist) - xlims[0]).argmin()
        right_idx = np.abs(np.array(proj_dist) - xlims[1]).argmin()
        proj_dist,shifted_mag = proj_dist[left_idx:right_idx],shifted_mag[left_idx:right_idx]
    zln = ax.plot(proj_dist,np.zeros(len(proj_dist)),linestyle=zlinestyle,color=zcolor,alpha=zalpha, clip_on=clip_on)
    sln = ax.plot(proj_dist,shifted_mag, clip_on=clip_on,**kwargs)

#    if not clip_on:
#        ax.set_xlim(xlims)
#        start_box = ax.transAxes.transform([0,0])[0]
#        end_box = ax.transAxes.transform([1,0])[0] - start_box
#        clip_patch = Rectangle([start_box,-1e7],end_box,1e9)
#        ax.lines[0].set_clip_path(clip_patch)
#        ax.lines[1].set_clip_path(clip_patch)

    if return_objects: return sln[0], zln[0]
    else: return min(proj_dist), max(proj_dist)

def remove_axis_lines_and_ticks(ax):
    [spline.set_visible(False) for spline in ax.spines.values()]
    ax.set_xticks([])
    ax.set_yticks([])

def plot_deskew_page(row, leave_plots_open=False, xlims=[-500,500], ylims=[-250,250], twf=0, **kwargs):
#    plt.rc('text', usetex=True)

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,8)
    ax0.set_anchor('W')
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%row['track_type'],rotation=0,fontsize=10)
    ax0.yaxis.set_label_coords(1.07,.45)
    if row['track_type']=="aero": layer_depth = 12.5
    else: layer_depth = 4.5
    min_syn_dis,max_syn_dis = plot_synthetic(row['sz_name'], row[['age_min','age_max']], ax0, layer_depth=layer_depth, color='k', linestyle='-', twf=twf)
    ax0.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
    ax0.patch.set_alpha(0.0)

    for j,phase_shift in enumerate(np.arange(row['phase_shift']-3*row['step'],row['phase_shift']+4*row['step'],row['step'])):
        if j>6: break #stop over plotting synthetic which can happen due to floats being used for phase shifts
        ax = fig.add_subplot(8,1,j+1, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        ax.set_ylabel(r"$\theta$=%.1f"%float(phase_shift),rotation=0,fontsize=14)
        ax.yaxis.set_label_coords(1.05,.45)
        min_proj_dis, max_proj_dis = plot_skewness_data(row,phase_shift,ax,picker=leave_plots_open)
        ax.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
        ax.patch.set_alpha(0.0)

    plot_chron_span_on_axes(row['sz_name'],fig.get_axes(),row[['age_min','age_max']])

#    ax0.set_xlim(max(min_proj_dis,min_syn_dis,-1000),min(max_proj_dis,max_syn_dis,1000))
    ax0.set_xlim(xlims)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    fig.suptitle("%s - PhaseShift=%.1f - Step=%d"%(row['comp_name'],float(row['phase_shift']),int(row['step'])),fontsize=16)
    plot_scale_bars(ax)

    out_file = os.path.join(row['results_dir'],'deskew_plots',"%s_%s_%s.png"%(str(row['comp_name']),str(round(row['phase_shift'],3)),str(row['step'])))

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open,msg="saving deskew file for %s to %s"%(row['comp_name'],out_file))

def plot_pole_perturbations(deskew_path,ellipse_path,**kwargs):

    deskew_df = filter_deskew_and_calc_aei(deskew_path)
    elipse_file = open(ellipse_path,"r")
    lon,lat,az,a,b = list(map(float,elipse_file.read().split()))

    for i,row in deskew_df.iterrows():
        run_in_parallel(plot_trial_pole_reduction_page,args=[row,lon,lat],kwargs=kwargs)

def plot_trial_pole_reduction_page(row,pole_lon, pole_lat, dis=5, spreading_rate_model_path=None, anomalous_skewness_model_path=None, xlims=[-500,500], ylims=[-200,200], leave_plots_open=False, twf=0, **kwargs):

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(6,1,6)
    ax0.set_anchor('W')
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%row['track_type'],rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(1.07,.45)
    if row['track_type']=="aero": layer_depth = 12.5
    else: layer_depth = 4.5
    min_syn_dis,max_syn_dis = plot_synthetic(row['sz_name'], row[['age_min','age_max']], ax0, layer_depth=layer_depth, color='k', linestyle='-', twf=twf)
    ax0.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
    ax0.patch.set_alpha(0.0)

    asf,srf,sz_list = get_asf_srf(spreading_rate_model_path,anomalous_skewness_model_path)
    azi = Geodesic.WGS84.Inverse(row['inter_lat'],row['inter_lon'],pole_lat,pole_lon)['azi1']

    for j in range(5):
        ax = fig.add_subplot(6,1,j+1, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        if j==0:
            geo_dict = Geodesic.WGS84.ArcDirect(pole_lat,pole_lon,dis,azi+90)
            tpole_lon,tpole_lat = geo_dict['lon2'],geo_dict['lat2']
            phase_shift = reduce_dsk_row_to_pole(row, tpole_lon, tpole_lat, asf, srf)
        elif j==1:
            geo_dict = Geodesic.WGS84.ArcDirect(pole_lat,pole_lon,dis,azi)
            tpole_lon,tpole_lat = geo_dict['lon2'],geo_dict['lat2']
            phase_shift = reduce_dsk_row_to_pole(row, tpole_lon, tpole_lat, asf, srf)
        elif j==2:
            tpole_lon,tpole_lat = pole_lon,pole_lat
            phase_shift = reduce_dsk_row_to_pole(row, tpole_lon, tpole_lat, asf, srf)
        elif j==3:
            geo_dict = Geodesic.WGS84.ArcDirect(pole_lat,pole_lon,dis,azi-180)
            tpole_lon,tpole_lat = geo_dict['lon2'],geo_dict['lat2']
            phase_shift = reduce_dsk_row_to_pole(row, tpole_lon, tpole_lat, asf, srf)
        elif j==4:
            geo_dict = Geodesic.WGS84.ArcDirect(pole_lat,pole_lon,dis,azi-90)
            tpole_lon,tpole_lat = geo_dict['lon2'],geo_dict['lat2']
            phase_shift = reduce_dsk_row_to_pole(row, tpole_lon, tpole_lat, asf, srf)
        else:
            raise ValueError("pole perturbation direction unknown")

        ax.set_ylabel(r"$\theta$=%.1f"%float(phase_shift),rotation=0,fontsize=14)
        ax.yaxis.set_label_coords(1.05,.45)
        ax.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
        ax.patch.set_alpha(0.0)


        min_proj_dis, max_proj_dis = plot_skewness_data(row,phase_shift,ax,picker=leave_plots_open)

    plot_chron_span_on_axes(row['sz_name'],fig.get_axes(),row[['age_min','age_max']])

    ax0.set_xlim(xlims)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    fig.suptitle("%s - PhaseShift=%.1f - Step=%d"%(row['comp_name'],float(row['phase_shift']),int(row['step'])),fontsize=16)
    plot_scale_bars(ax)

    check_dir(os.path.join(row['results_dir'],'pole_perturbation_plots'))

    out_file = os.path.join(row['results_dir'],'pole_perturbation_plots',"%s_%s_%s.png"%(str(row['comp_name']),str(round(row['phase_shift'],3)),str(row['step'])))

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open,msg="saving deskew file for %s to %s"%(row['comp_name'],out_file))

def save_show_plots(fig,out_file,leave_plots_open=True,msg=None):
    if isinstance(msg,str): print(msg)
    if leave_plots_open:
        plt.show()
    else:
        fig.savefig(out_file,transparent=False,facecolor="white")
        plt.close(fig)

def plot_skewnesses(deskew_path,**kwargs):
    deskew_df = open_deskew_file(deskew_path)

    check_dir(os.path.join(deskew_df.iloc[0]['results_dir'],'deskew_plots'))

    for i,row in deskew_df.iterrows():
        print("starting process for %s"%row['comp_name'])
        run_in_parallel(plot_deskew_page,args=[row],kwargs=kwargs)
#        plot_deskew_page(row,**kwargs)

def plot_ridge_loc(deskew_row,ridge_loc_func,ax,**kwargs):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = open_mag_file(data_file_path)
    track_lon_lats = [[convert_to_0_360(lon),lat] for lon,lat in zip(data_df['lon'].tolist(),data_df['lat'].tolist())]

    ridge_lon,ridge_lat = ridge_loc_func(deskew_row['sz_name'],track_lon_lats)

    if ridge_lon==None or ridge_lat==None: print("ridge intercept values given are None"); return

    ridge_dict = Geodesic.WGS84.Inverse(float(deskew_row['inter_lat']),convert_to_0_360(deskew_row['inter_lon']),ridge_lat,ridge_lon)
    ridge_dis = (ridge_dict['s12']*np.sin(np.deg2rad(float(deskew_row['strike'])-ridge_dict['azi2'])))/1000

    ax.axvline(ridge_dis,**kwargs)

def plot_fz_loc(deskew_row,fz_loc_path,ax,**kwargs):
    fzdf = pd.read_csv(fz_loc_path,sep='\t',index_col=0)

    try:
        fz_inters = fzdf['inters'][deskew_row['comp_name']]
        fz_inters = list(map(lambda x: list(map(lambda y: float(y.strip(" '")), x.strip(", ] ['").split(','))), fz_inters.split(']')[:-2]))
    except KeyError: return

    for fz_lon,fz_lat in fz_inters:
        fzi_dict = Geodesic.WGS84.Inverse(float(deskew_row['inter_lat']),convert_to_0_360(deskew_row['inter_lon']),fz_lat,fz_lon)
        fz_dis = (fzi_dict['s12']*np.sin(np.deg2rad(float(deskew_row['strike'])-fzi_dict['azi2'])))/1000
        ax.axvline(fz_dis,**kwargs)

def plot_best_skewness_page(rows,results_dir,page_num,leave_plots_open=False,ridge_loc_func=None,fz_loc_path=None, xlims=[-500,500], ylims=[-250,250], clip_on = False, twf=0, layer_mag=1000, **kwargs):
#    plt.rc('text', usetex=True)
#    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

    fig = plt.figure(figsize=(12, 9), facecolor='white')

#    ax0 = fig.add_subplot(len(rows)+2,1,1)
#    remove_axis_lines_and_ticks(ax0)
#    ax0.set_ylabel(r"synthetic (%s)"%'ship',rotation=0,fontsize=10)
#    ax0.yaxis.set_label_coords(-.075,.45)
#    ax0.format_coord = format_coord
#    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax0, layer_depth=4.5, color='k', linestyle='-', clip_on=clip_on, xlims=xlims, twf=twf, layer_mag=layer_mag)
##    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax0, layer_depth=12.5, color='k', linestyle='-', clip_on=clip_on, xlims=xlims)
#    ax0.set_ylim(ylims) #MODIFY THIS TO CHANGE Y AXIS
#    ax0.patch.set_alpha(0.0)
#    ylim = ax0.get_ylim()

    for j,(i,row) in enumerate(rows.iterrows()):
#        ax = fig.add_subplot(len(rows)+2,1,j+2, sharex=ax0)
        if j==0:
            ax0 = fig.add_subplot(len(rows),1,j+1)
            ax = ax0
        else:
            ax = fig.add_subplot(len(rows),1,j+1, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        min_proj_dis, max_proj_dis = plot_skewness_data(row,float(row['phase_shift']),ax,picker=True, clip_on=clip_on, xlims=xlims)

        if row["track_type"]=="aero":
            try: old_data_df = utl.open_mag_file(os.path.join(row["data_dir"],row["comp_name"].replace(".Ed","").replace(".Vd","").replace(".Hd","")))
            except FileNotFoundError:
                if sys.platform=="win32": sep = "\\"
                else: sep = "/"
                cut_number = row["comp_name"].split(".")[1]
                old_data_df = utl.open_mag_file(os.path.join(*row["data_dir"].strip(sep).split(sep)[:-1],row["comp_name"].replace(".%s.Ed.lp"%cut_number,".DAT").replace(".%s.Vd.lp"%cut_number,".DAT").replace(".%s.Hd.lp"%cut_number,".DAT")))
            tmp_projected_distances = utl.calc_projected_distance(row['inter_lon'],row['inter_lat'],old_data_df['lon'].tolist(),old_data_df['lat'].tolist(),(180+row['strike'])%360)
            anomaly_middle = np.argwhere(np.diff(np.sign(tmp_projected_distances["dist"])))[0] #this is not to only get the first 0 there is only 1 it's because of a wrapper sequence
            layer_depth = (4.5 + 0.0003048*(old_data_df["alt"][anomaly_middle]))
            if isinstance(layer_depth,pd.Series): layer_depth = layer_depth.iloc[0]
        else: layer_depth = 4.5 #default approx depth to layer 2A in deep Pacific
        min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax, layer_depth=layer_depth, color='r', linestyle='-', linewidth=2, alpha=.5, clip_on=clip_on, xlims=xlims, twf=twf, layer_mag=layer_mag, zorder=-1)

        if ridge_loc_func!=None:
            plot_ridge_loc(row,ridge_loc_func,ax,color='r',linestyle='-',alpha=1)
        if fz_loc_path!=None:
            plot_fz_loc(row,fz_loc_path,ax,color='g',linestyle='-',alpha=1)

        ax.annotate(r"%s"%row['comp_name']+"\n"+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(row['inter_lat']),convert_to_0_360(row['inter_lon'])),xy=(-.15,.45),xycoords="axes fraction",fontsize=10)
        ax.set_ylabel(r"$\theta$=%.1f"%float(row['phase_shift'])+"\n"+r"$e_a$=%.1f"%float(row['aei']),rotation=0,fontsize=10)
        ax.yaxis.set_label_coords(1.05,.45)
        ax.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
        ax.patch.set_alpha(0.0)
        ax.format_coord = format_coord
#        for l in ax.lines:
#            import pdb; pdb.set_trace()
#            print(ax.transData.inverted().transform(l.get_clip_path().get_xy()),ax.transData.inverted().transform([l.get_clip_path().get_width(),l.get_clip_path().get_height()]))
#        print("----")

    for j in range(j+1,len(rows)):
        ax = fig.add_subplot(len(rows)+2,1,j+2, sharex=ax0)
        ax.set_anchor('W')
        remove_axis_lines_and_ticks(ax)
        ax.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
        ax.patch.set_alpha(0.0)
        ax.format_coord = format_coord

#    ax = fig.add_subplot(len(rows)+2,1,len(rows)+2, sharex=ax0)
#    ax.set_anchor('W')
#    remove_axis_lines_and_ticks(ax)
#    ax.set_ylabel(r"synthetic (%s)"%'aero',rotation=0,fontsize=10)
#    ax.yaxis.set_label_coords(-.075,.45)
#    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax, layer_depth=12.5, color='k', linestyle='-', clip_on=clip_on, xlims=xlims, twf=twf, layer_mag=layer_mag)
#    ax.set_ylim(ylims) #insure that all of the plots have the same zoom level in the y direction
#    ax.patch.set_alpha(0.0)
#    ax.format_coord = format_coord

    plot_chron_span_on_axes(rows['sz_name'].iloc[0],fig.get_axes(),rows[['age_min','age_max']].iloc[0])

    ax0.set_xlim(xlims)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    if rows.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
        title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows['sz_name'].iloc[0],float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    else:
        title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    fig.suptitle(title,fontsize=16)
    plot_scale_bars(ax)

    out_file = os.path.join(results_dir,"page_%d.pdf"%page_num)

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open)

def plot_best_skewnesses(deskew_df, best_skews_subdir="best_skews", num_profiles_per_page = 6, **kwargs):
    results_dir = os.path.join(deskew_df['results_dir'].iloc[0],best_skews_subdir)
    check_dir(results_dir)

    if "ridge_loc_path" in kwargs.keys() and kwargs["ridge_loc_path"]!=None: ridge_loc_func = read_and_fit_ridge_data(kwargs["ridge_loc_path"])

    prev_i,page_num = 0,0
    for i in list(range(num_profiles_per_page,len(deskew_df.index),num_profiles_per_page))+[-1]:
        if i==-1: rows = deskew_df.iloc[prev_i:]
        else: rows = deskew_df.iloc[prev_i:i]
        run_in_parallel(plot_best_skewness_page,args=[rows,results_dir,page_num],kwargs=kwargs)
#        plot_best_skewness_page(rows,results_dir,page_num,**kwargs)
        prev_i,page_num = i,page_num+1

def overlay_skewness_page(rows1,rows2,results_dir,page_num,leave_plots_open=False,pole_name1='pole 1', pole_name2='pole 2', fz_loc_path=None, twf=0):
#    plt.rc('text', usetex=True)
#    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,1)
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%'ship',rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(-.1,.45)
    ax0.format_coord = format_coord
    min_syn_dis,max_syn_dis = plot_synthetic(rows1['sz_name'].iloc[0], rows1[['age_min','age_max']].iloc[0], ax0, layer_depth=4.5, color='k', linestyle='-', twf=twf)
    ylim = ax0.set_ylim(-200,200)

    for j,iterrow in enumerate(zip(rows1.iterrows(),rows2.iterrows())):
        iterrow1,iterrow2 = iterrow
        i,row1 = iterrow1
        k,row2 = iterrow2
        ax = fig.add_subplot(8,1,j+2, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        min_proj_dis, max_proj_dis = plot_skewness_data(row1,float(row1['phase_shift']),ax,color='k',linestyle='-',picker=True,alpha=.5)
        min_proj_dis, max_proj_dis = plot_skewness_data(row2,float(row2['phase_shift']),ax,color='r',linestyle='-',picker=True,alpha=.5)

        if fz_loc_path!=None:
            plot_fz_loc(row1,fz_loc_path,ax,color='g',linestyle='-',alpha=1)

        ax.annotate(r"%s"%row1['comp_name']+"\n"+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(row1['inter_lat']),convert_to_0_360(row1['inter_lon'])),xy=(-.15,.45),xycoords="axes fraction")
        ax.set_ylabel(r"$\theta$=%.1f"%float(row1['phase_shift'])+"\n"+r"$e_a$=%.1f"%float(row1['aei'])+'\n'+r"$\theta_2$=%.1f"%float(row2['phase_shift'])+"\n"+r"$e_{a2}$=%.1f"%float(row2['aei'])+'\n-----------',rotation=0,fontsize=10)
        ax.yaxis.set_label_coords(1.05,0.0)
        ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction
        ax.format_coord = format_coord

    ax = fig.add_subplot(8,1,8, sharex=ax0)
    ax.set_anchor('W')
    remove_axis_lines_and_ticks(ax)
    ax.set_ylabel("synthetic (%s)"%'aero',rotation=0,fontsize=14)
    ax.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows1['sz_name'].iloc[0], rows1[['age_min','age_max']].iloc[0], ax, layer_depth=12.5, color='k', linestyle='-', twf=twf)
    ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction
    ax.format_coord = format_coord

    plot_chron_span_on_axes(rows1['sz_name'].iloc[0],fig.get_axes(),rows1[['age_min','age_max']].iloc[0])

    ax0.set_xlim(-400,600)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes

    if rows1.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
        title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows1['sz_name'].iloc[0],float(rows1['inter_lat'].iloc[0]),float(rows1['inter_lat'].iloc[-1]))
    else:
        title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows1['inter_lat'].iloc[0]),float(rows1['inter_lat'].iloc[-1]))
    fig.suptitle(title,fontsize=16)
    plot_scale_bars(ax)

    poles_handles=[]
    for pole_name,color,style in [[pole_name1,'k','-'],[pole_name2,'r','-']]:
        poles_handle = mlines.Line2D([], [], color=color, linestyle=style, label=pole_name)
        poles_handles.append(poles_handle)

    plt.legend(handles=poles_handles,loc=3,bbox_to_anchor=(.0, .0),bbox_transform=plt.gcf().transFigure, frameon=False)

    out_file = os.path.join(results_dir,"page_%d.png"%page_num)

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open)

def overlay_best_skewnesses(deskew_path, deskew_path2, leave_plots_open=False, best_skews_subdir="best_skews", pole_name1='pole 1', pole_name2='pole 2', fz_loc_path=None, num_profiles_per_page = 6):
#    dt_path,age_min,age_max,results_dir = create_matlab_datatable(deskew_path)

    deskew_df = filter_deskew_and_calc_aei(deskew_path)
    deskew_df2 = filter_deskew_and_calc_aei(deskew_path2)

    results_dir = os.path.join(deskew_df['results_dir'].iloc[0],best_skews_subdir)
    check_dir(results_dir)

    prev_i,page_num = 0,0
    for i in range(num_profiles_per_page,len(deskew_df.index),num_profiles_per_page):
        rows = deskew_df.iloc[prev_i:i]
        rows2 = deskew_df2.iloc[prev_i:i]
        page_num = i/num_profiles_per_page
        run_in_parallel(overlay_skewness_page,args=[rows,rows2,results_dir,page_num],kwargs={'leave_plots_open':leave_plots_open,'pole_name1':pole_name1,'pole_name2':pole_name2,'fz_loc_path':fz_loc_path})
        prev_i = i
    rows = deskew_df.iloc[prev_i:len(deskew_df.index)]
    rows2 = deskew_df2.iloc[prev_i:len(deskew_df2.index)]
    page_num += 1
    run_in_parallel(overlay_skewness_page,args=[rows,rows2,results_dir,page_num],kwargs={'leave_plots_open':leave_plots_open,'pole_name1':pole_name1,'pole_name2':pole_name2,'fz_loc_path':fz_loc_path})

def overlay_skewness_by_spreading_zone(deskew_path,deskew_path2,leave_plots_open=False,pole_name1='pole 1', pole_name2='pole 2', fz_loc_path=None):
    deskew_df = open_deskew_file(deskew_path)
    deskew_df2 = pd.read_csv(deskew_path2,sep='\t')

    tmp_deskew_path = '.tmp_deskew_file'
    tmp_deskew_path2 = '.tmp2_deskew_file'
    for sz in deskew_df['sz_name'].drop_duplicates():
        sz_deskew_df = deskew_df[deskew_df['sz_name']==sz].copy()
        sz_deskew_df2 = deskew_df2[deskew_df2['sz_name']==sz].copy()
        sz_deskew_df.sort_values(by="comp_name",inplace=True)
        sz_deskew_df2.sort_values(by="comp_name",inplace=True)
        sz_deskew_df.sort_values(by="inter_lat",ascending=False,inplace=True)
        sz_deskew_df2.sort_values(by="inter_lat",ascending=False,inplace=True)
        sz_deskew_df.to_csv(tmp_deskew_path,sep='\t',index=False)
        sz_deskew_df2.to_csv(tmp_deskew_path2,sep='\t',index=False)
        overlay_best_skewnesses(tmp_deskew_path,tmp_deskew_path2,leave_plots_open=leave_plots_open,best_skews_subdir="overlayed_skewness_by_spreading_zones/%s"%str(sz), pole_name1=pole_name1, pole_name2=pole_name2, fz_loc_path=fz_loc_path)

    try: os.remove(tmp_deskew_path); os.remove(tmp_deskew_path2); os.remove("%s.dt.csv"%os.path.basename(tmp_deskew_path).split('.')[0])
    except OSError: print("trouble removing temporary deskew and dt.csv files used for this function, check the code")

def plot_skewness_by_spreading_zone(deskew_path,leave_plots_open=False,ridge_loc_path=None, fz_loc_path=None, synth_config_path=None, **kwargs):
    deskew_df = filter_deskew_and_calc_aei(deskew_path)

    for sz in deskew_df['sz_name'].drop_duplicates():
        sz_deskew_df = deskew_df[deskew_df['sz_name']==sz].copy()
        sz_deskew_df.sort_values(by="comp_name",inplace=True)
        sz_deskew_df.sort_values(by="inter_lat",ascending=False,inplace=True)
        plot_best_skewnesses(sz_deskew_df,leave_plots_open=leave_plots_open,ridge_loc_path=ridge_loc_path,best_skews_subdir="skewness_by_spreading_zones/%s"%str(sz), fz_loc_path=fz_loc_path, **kwargs)

def plot_spreading_rate_picks_page(rows, spreading_rate_picks, results_dir, page_num, leave_plots_open=False, twf=0):

    check_dir(results_dir)

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,1)
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%'ship',rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax0, layer_depth=4.5, color='k', linestyle='-', plot_anom_spans=True, twf=twf)
    ylim = ax0.get_ylim()

    for j,iterrow in enumerate(rows.iterrows()):
        i,row = iterrow
        ax = fig.add_subplot(8,1,j+2, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        if row['track_type'] == 'aero':
            row['comp_name'] += '.Ed.lp'
        elif row['track_type'] == 'ship':
            row['comp_name'] += '.lp'
        else: print("could not determine data type (ship or aero) of %s, skipping"%row['comp_name']); continue

        min_proj_dis, max_proj_dis = plot_skewness_data(row,float(row['phase_shift']),ax,picker=True)

        ax.annotate(r"%s"%row['comp_name'].rstrip('.Ed .Vd .lp')+"\n"+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(row['inter_lat']),convert_to_0_360(row['inter_lon'])),xy=(-.15,.45),xycoords="axes fraction")
        ax.set_ylabel(r"$\theta$=%.1f"%float(row['phase_shift']),rotation=0,fontsize=14)
        ax.yaxis.set_label_coords(1.05,.45)
        ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

        for anomoly,pick in spreading_rate_picks[row['comp_name'].rstrip('.Ed .Vd .lp')].iteritems():
            if np.isnan(pick): continue
            ax.axvline(pick,linestyle='--',color='red',alpha=.7)
            if 'n' in anomoly.split('.')[-1]: ylevel=ax.get_ylim()[0]-ax.get_ylim()[0]/5
            else: ylevel=ax.get_ylim()[0]
            ax.annotate(anomoly, xy=(pick,ylevel), xycoords='data', ha='left', va='bottom', fontsize=8, clip_on=False)

    ax = fig.add_subplot(8,1,8, sharex=ax0)
    ax.set_anchor('W')
    remove_axis_lines_and_ticks(ax)
    ax.set_ylabel("synthetic (%s)"%'aero',rotation=0,fontsize=14)
    ax.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0], rows[['age_min','age_max']].iloc[0], ax, layer_depth=12.5, color='k', linestyle='-', plot_anom_spans=True, twf=twf)
    ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

    ax0.set_xlim(-1000,500)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    if rows.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
        title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows['sz_name'].iloc[0],float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    else:
        title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    fig.suptitle(title,fontsize=16)
    plot_scale_bars(ax)

    out_file = os.path.join(results_dir,"page_%d.png"%page_num)

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open,msg="Saving as %s"%out_file)

def plot_spreading_rate_picks(deskew_path,spreading_rate_picks_path,leave_plots_open=False):
    deskew_df = open_deskew_file(deskew_path)
    spreading_rate_picks = pd.read_csv(spreading_rate_picks_path,sep='\t',header=0,index_col=0)

    deskew_df['comp_name'] = deskew_df['comp_name'].apply(lambda x: str(x).rstrip('.Ed .Vd .lp'))
    picked_profiles = deskew_df[deskew_df['comp_name'].isin(spreading_rate_picks.columns)].drop_duplicates(subset='comp_name')

    results_dir = os.path.join(deskew_df['results_dir'].iloc[0],'spreading_rate_picks')

    num_profiles_per_page = 6
    prev_i,page_num = 0,0
    for i in range(num_profiles_per_page,len(picked_profiles.index),num_profiles_per_page):
        rows = picked_profiles.iloc[prev_i:i]
        page_num = i/num_profiles_per_page
        plot_spreading_rate_picks_page(rows, spreading_rate_picks, results_dir, page_num, leave_plots_open=leave_plots_open)
        prev_i = i
    rows = picked_profiles.iloc[prev_i:len(picked_profiles.index)]
    page_num += 1
    plot_spreading_rate_picks_page(rows, spreading_rate_picks, results_dir, page_num, leave_plots_open=leave_plots_open)

def plot_spreading_rate_results(sr_path, xlims=[-500,500], ylims=[-200,200], median_or_mean='mean', title="", leave_plots_open=False):
    sr_df = pd.read_csv(sr_path,sep='\t',header=0,index_col=0)

    mean_median_values = [val for val in sr_df[median_or_mean].tolist() for _ in (0, 1)]
    std_values = [val for val in sr_df['std'].tolist() for _ in (0, 1)]
    n_values = [val for val in sr_df['n'].tolist() for _ in (0, 1)]
    nested_ages = [(min_a,max_a) for min_a,max_a in zip(sr_df['age_min'].tolist(),sr_df['age_max'].tolist())]
    ages = [item for sublist in nested_ages for item in sublist]

    confidence_ub = np.array(mean_median_values)+(np.array(std_values)*2.0)/np.sqrt(np.array(n_values))
    confidence_lb = np.array(mean_median_values)+(np.array(std_values)*-2.0)/np.sqrt(np.array(n_values))

    fig = plt.figure(figsize=(16, 9))
    ax = plt.gca()
    ax.invert_xaxis()

    plt.plot(ages,mean_median_values,color='black',linewidth=3)
    plt.plot(ages,confidence_ub,color='black',linestyle='--',linewidth=1)
    plt.plot(ages,confidence_lb,color='black',linestyle='--',linewidth=1)

    # Fix x-axis limits if specified
    plt.xlim( xlims )
    plt.ylim( ylims )

    i = 0
    for mean_median,std,age,n in zip(sr_df[median_or_mean].tolist(),sr_df['std'].tolist(),nested_ages,sr_df['n'].tolist()):
        anom = sr_df.index.where(sr_df['std'] == std)[i]
        ax.add_patch(Rectangle((age[0], mean_median - (2.0*std)/np.sqrt(n)), age[1]-age[0],  (4.0*std)/np.sqrt(n),color='grey',alpha=.5,zorder=0,linewidth=0))
        ax.annotate('%s'%anom,xy=(sum(age)/len(age),mean_median + (2.0*std)/np.sqrt(n)),va='bottom',ha='center',fontsize=8)
        ax.annotate('n=%d'%n,xy=(sum(age)/len(age),mean_median - (2.0*std)/np.sqrt(n) - 1),va='top',ha='center',fontsize=8)

        i += 1

    plt.xlabel("Age (Myr)")
    plt.ylabel("Spreading Half Rate (km/Myr)")
    plt.title(title)

    out_file = os.path.join(os.path.dirname(sr_path),os.path.basename(sr_path).split('.')[0]+"_spreading_rate_%s.png"%(median_or_mean))

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open)

def plot_isochron_picks(deskew_path,spreading_rate_picks_path,leave_plots_open=False):
    iso_df,average_lon,average_lat = get_lon_lat_from_plot_picks_and_deskew_file(deskew_path,spreading_rate_picks_path)
    if iso_df.empty: print("problem with input data, aborting"); return
    average_lon = convert_to_0_360(average_lon)

    # Read data files
    deskew_df = open_deskew_file(deskew_path)
    results_dir = deskew_df['results_dir'].iloc[0]

    # Initialize figure
    fig = plt.figure(figsize=(16,9), dpi=200)

    lon_spread,lat_spread = 20,15

    llcrnrlon=round(average_lon/10)*10-lon_spread
    llcrnrlat=round(average_lat/10)*10-lat_spread
    urcrnrlon=round(average_lon/10)*10+lon_spread
    urcrnrlat=round(average_lat/10)*10+lat_spread

    m = create_basic_map(projection='merc', llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    m.drawmeridians(np.arange(llcrnrlon,urcrnrlon-5,5),labels=[0,0,0,1],color=[0,0,0,.1],zorder=0)
    m.drawparallels(np.arange(llcrnrlat,urcrnrlat+5,5),labels=[1,0,0,0],color=[0,0,0,.1],zorder=0)

    sort_lonlat_cmp = lambda x,y=None: x[1]-y[1] if not isinstance(y,type(None)) else -1
    sort_lonlat_key = cmp_to_key(sort_lonlat_cmp)

    chrons_info = []
    for anom,row in iso_df.iterrows():
        row = row[row.notnull()]
        if len(row)<3: print("cannot plot isochron for %s as it only has less than 3 data points so it will be skipped"%anom); continue
        lonlats = [[convert_to_0_360(latlon['lon']),float(latlon['lat'])] for latlon in row.tolist()]
        lonlats.sort(key=sort_lonlat_key)
        lonlats = np.array(lonlats)
        iso_line, = m.plot(lonlats[:,0], lonlats[:,1],  label=anom, zorder=3, marker="x")
        if anom.strip('C r n') not in list(map(lambda x: x[0], chrons_info)): chrons_info.append([anom.strip('C r n'), iso_line.get_color()])

    plot_chron_info(chrons_info,m,coord_0_360=True,zorder=2,alpha=.4,linestyle='--',dashes=(5,1))

    legend = plt.legend()
    frame = legend.get_frame()
    frame.set_alpha(.7)

    out_file = os.path.join(results_dir,'iso_chron_plots',os.path.basename(spreading_rate_picks_path).split('.')[0]+'.png')

    save_show_plots(fig,out_file,leave_plots_open=leave_plots_open,msg="Saving to %s"%out_file)


