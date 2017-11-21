import os
import subprocess
import pandas as pd
import numpy as np
import pmagpy.ipmag as ipmag
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import OrderedDict
from matplotlib.patches import Rectangle
from multiprocessing import Process
from geographiclib.geodesic import Geodesic
from .skewness import *
from .plot_geographic import *
from .utilities import *

def plot_all_lunes_seperate(deskew_path):

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
            fig = plt.figure(figsize=(16,9), dpi=80)

            # Create figure
            gcm = create_basic_map(projection='npstere', lat_0=90, boundinglat=60)
            gcm.drawparallels(np.arange(0,90,10),labels=[0,0,0,0])
            gcm.drawmeridians(np.arange(0,360,10),labels=[1,1,0,1])

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
        gcm.plot(gc_lon, gc_lat, color='black',latlon=True,linestyle=linestyle,label=label)

        if row["track_type"]=='ship' or last_comp_dir in row["comp_name"]:
            plt.title(track_name)
            plt.legend(loc='best')

            plot_path = os.path.join(lunes_dir,track_name+".png")

            print("saving to %s"%plot_path)
            fig.savefig(plot_path)
            plt.close(fig)

def plot_pole(ellipse_path,m=None,color='cyan',marker='o',s=30,zorder=3,alpha=.5,label=None,**kwargs):

    elipse_file = open(ellipse_path,"r")
    lon,lat,az,a,b = list(map(float,elipse_file.read().split()))

    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=80)
        ax = fig.add_subplot(111)
        # Create map
        m = create_basic_map(projection='npstere', lat_0=90, boundinglat=60, ax=ax)
        m.drawparallels(np.arange(0,90,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,1,1])

    m.scatter(lon, lat, facecolor=color, edgecolor='black', marker=marker, s=s, latlon=True, zorder=zorder, label=label)
    m.scatter(0, 90, facecolor='black', edgecolor='black', marker='+', s=s, latlon=True, zorder=1, label=label)
    ipmag.ellipse(m, lon, lat, (a*111.11)/2, (b*111.11)/2, az, n=360, filled=True, facecolor=color, edgecolor='black', zorder=zorder-1,alpha=alpha)

    return m,lon,lat,az,a,b

def plot_pole_track(ellipse_paths, m=None, **kwargs):

    if m==None:
        # Create figure
        fig = plt.figure(figsize=(16,9), dpi=80)
        ax = fig.add_subplot(111)
        # Create map
        m = create_basic_map(projection='npstere', lat_0=90, boundinglat=60, ax=ax)
        m.drawparallels(np.arange(0,90,10),labels=[0,0,0,0])
        m.drawmeridians(np.arange(0,360,10),labels=[1,1,1,1])

    lats,lons = [],[]
    if 'label' in kwargs.keys(): label = kwargs.pop('label')
    for ellipse_path in ellipse_paths:
        m,lon,lat,az,a,b = plot_pole(ellipse_path,m=m,**kwargs)
        lons.append(lon);lats.append(lat)

    m.plot(lons,lats,latlon=True,label=label,**kwargs)

    return m

def plot_pole_with_lunes(deskew_path,ellipse_path):

    comps = filter_deskew_and_calc_aei(deskew_path)

    # Create figure
    fig = plt.figure(figsize=(16,9), dpi=80)
    ax = fig.add_subplot(111)

    # Create map
    gcm = create_basic_map(projection='npstere', lat_0=90, boundinglat=60, ax=ax)
    gcm.drawparallels(np.arange(0,90,10),labels=[0,0,0,0])
    gcm.drawmeridians(np.arange(0,360,10),labels=[1,1,1,1])

    plot_lunes(comps,gcm)

    plot_pole(ellipse_path,m=gcm)

    out_path = os.path.join(comps['results_dir'].iloc[0],"poles.png")
    check_dir(os.path.dirname(out_path))
    print("saving to %s"%out_path)
    fig.savefig(out_path)
    plt.close(fig)

def plot_lunes_and_save(deskew_path):

    comps = filter_deskew_and_calc_aei(deskew_path)

    # Create figure
    fig = plt.figure(figsize=(16,9), dpi=80)

    # Create figure
    gcm = create_basic_map(projection='npstere', lat_0=90, boundinglat=60)
    gcm.drawparallels(np.arange(0,90,10),labels=[0,0,0,0])
    gcm.drawmeridians(np.arange(0,360,10),labels=[1,1,1,1])

    plot_lunes(comps,gcm)

    out_path = os.path.join(comps['results_dir'].iloc[0],"lunes.png")
    print("saving to %s"%out_path)
    fig.savefig(out_path)
    plt.close(fig)


def plot_lunes(comps,gcm):

    # For every crossing...
    for i,row in comps.iterrows():
        if row["comp_name"].startswith('#'): continue

        if row["track_type"]=='aero':
            if "V" in row["comp_name"]: linestyle="--"
            else: linestyle="-"
            linewidth=1
        if row["track_type"]=='ship':
            linestyle,linewidth=":",1

        strike = convert_to_0_360(row['strike'])
        inc_mid = np.degrees(np.arctan(np.tan(np.deg2rad(float(row["aei"])))))
        clt_mid = 90-np.degrees(np.arctan(np.tan(np.deg2rad(inc_mid))*0.5))
        lat_mid = Geodesic.WGS84.ArcDirect(float(row["inter_lat"]),float(row["inter_lon"]), strike-90, clt_mid)['lat2']
        # Find array of great semicircle azimuths (degrees)
        if lat_mid >= 0:
            azi = np.linspace(strike-(180),strike,100)
        else:
            azi = np.linspace(strike,strike+180,100)

        # For first bounding azimuth...
        # Find inclination (degrees)
        inc = np.degrees(np.arctan(np.tan(np.deg2rad(float(row["aei"])))*np.sin(np.deg2rad(azi+180-strike))))
        # Find paleocolatitude (degrees)
        clt = 90-np.degrees(np.arctan(np.tan(np.deg2rad(inc))*0.5))
        # Find great circle points
        gc_points_and_azis = [Geodesic.WGS84.ArcDirect(float(row["inter_lat"]),float(row["inter_lon"]), azi1, clt1) for clt1,azi1 in zip(clt,azi)]
        gc_lon = [gcd["lon2"] for gcd in gc_points_and_azis]
        gc_lat = [gcd["lat2"] for gcd in gc_points_and_azis]

        # Draw great circle
        gcm.plot(gc_lon, gc_lat, color=(float(row["r"]),float(row["g"]),float(row["b"])),linestyle=linestyle,linewidth=linewidth,latlon=True)

    comps["inter_lat"] = comps["inter_lat"].apply(float)
    tmp_comps = comps.sort_values(by="inter_lat",ascending=False)

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

    track_type_legend = plt.legend(handles=track_type_handles,loc=2)
    tt_frame = track_type_legend.get_frame()
    tt_frame.set_alpha(.7)
    sz_legend = plt.legend(handles=sz_handles,bbox_to_anchor=(.95, .91),bbox_transform=plt.gcf().transFigure, frameon=False)
    sz_frame = sz_legend.get_frame()
    sz_frame.set_alpha(.7)

    plt.gca().add_artist(track_type_legend)

    return gcm

def plot_chron_span_on_axes(sz,axes):
    # I know I made it save two files, but that was before I remembered it should be the same. :p
    try: dis_span = np.loadtxt(os.path.join('SynthData','dis_span_%s.txt'%sz),delimiter=',')
    except (IOError,OSError) as e:
        try:

            default_name = ''
            if os.path.isfile(os.path.join('SynthData','dis_span_Default.txt')): default_name = 'Default'
            elif os.path.isfile(os.path.join('SynthData','dis_span_default.txt')): default_name = 'default'
            else: sz = 'default'; raise IOError

            dis_span = np.loadtxt(os.path.join('SynthData','dis_span_%s.txt'%default_name),delimiter=',')

        except (IOError,OSError) as e:
            raise IOError("No synthetic found for %s, please generate a new synthetic which contains a spreading rate model for this spreading zone"%sz)

    for axis in axes:
        axis.axvspan(*dis_span, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)

def plot_synthetic(sz,ship_or_aero,ax,**kwargs):
    try:
        if ship_or_aero=='aero':
            dis_syn = np.loadtxt(os.path.join('SynthData','dis_syn_%s_aero.txt'%sz),delimiter=',')
            mag_syn = np.loadtxt(os.path.join('SynthData','mag_syn_%s_aero.txt'%sz),delimiter=',')
        elif ship_or_aero=='ship':
            dis_syn = np.loadtxt(os.path.join('SynthData','dis_syn_%s_ship.txt'%sz),delimiter=',')
            mag_syn = np.loadtxt(os.path.join('SynthData','mag_syn_%s_ship.txt'%sz),delimiter=',')
        else: raise ValueError("plot synthetic needs to know if you want the aeromag or shipmag synthetic was given %s not ship or aero"%ship_or_aero)
    except (IOError,OSError) as e:
        try:

            default_name = ''
            if os.path.isfile(os.path.join('SynthData','dis_span_Default.txt')): default_name = 'Default'
            elif os.path.isfile(os.path.join('SynthData','dis_span_default.txt')): default_name = 'default'
            else: sz = 'default'; raise IOError

            if ship_or_aero=='aero':
                dis_syn = np.loadtxt(os.path.join('SynthData','dis_syn_%s_aero.txt'%default_name),delimiter=',')
                mag_syn = np.loadtxt(os.path.join('SynthData','mag_syn_%s_aero.txt'%default_name),delimiter=',')
            elif ship_or_aero=='ship':
                dis_syn = np.loadtxt(os.path.join('SynthData','dis_syn_%s_ship.txt'%default_name),delimiter=',')
                mag_syn = np.loadtxt(os.path.join('SynthData','mag_syn_%s_ship.txt'%default_name),delimiter=',')
            else: raise ValueError("plot synthetic needs to know if you want the aeromag or shipmag synthetic was given %s not ship or aero"%ship_or_aero)

        except (IOError,OSError) as e:
            raise IOError("No synthetic found for %s, please generate a new synthetic which contains a spreading rate model for this spreading zone"%sz)

    if 'plot_anom_spans' in kwargs and kwargs.pop('plot_anom_spans'):
        if not os.path.isfile(os.path.join('SynthData','anom_spans_%s.txt'%sz)):
            print("%s does not exist and is required in order to plot anomoly spans on the synthetics, please regenerate the synthetic with an age model file to solve this issue."%os.path.abspath(os.path.join('SynthData','anom_spans.txt')))
        else:
            anom_spans = pd.read_csv(os.path.join('SynthData','anom_spans_%s.txt'%sz),sep=',',index_col=0,header=None)
            for anomoly,older_bound in anom_spans.iterrows():
                if isinstance(older_bound,pd.Series): older_bound=older_bound.iloc[0]
                if np.isnan(older_bound): continue
                ax.axvline(older_bound,linestyle='--',color='blue',alpha=.7)
                if 'n' in anomoly.split('.')[-1]: ylevel=-380+80
                else: ylevel=-380
                ax.annotate(anomoly, xy=(older_bound,ylevel), xycoords='data', ha='left', va='bottom', fontsize=8, clip_on=False)

    ax.plot(dis_syn,np.zeros(len(dis_syn)),'k--')
    ax.plot(dis_syn,mag_syn,**kwargs)

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

def plot_skewness_data(deskew_row, phase_shift, ax, **kwargs):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","dec_year","mag","lat","lon"],delim_whitespace=True)

    projected_distances = calc_projected_distance(deskew_row['inter_lon'],deskew_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),deskew_row['strike'])

    shifted_mag = phase_shift_data(data_df['mag'].tolist(),phase_shift)

    proj_dist = projected_distances['dist'].tolist()
    ax.plot(proj_dist,np.zeros(len(proj_dist)),'k--')
    ax.plot(proj_dist,shifted_mag,**kwargs)

    return min(proj_dist), max(proj_dist)

def remove_axis_lines_and_ticks(ax):
    [spline.set_visible(False) for spline in ax.spines.values()]
    ax.set_xticks([])
    ax.set_yticks([])

def plot_deskew_file_skewnesses(row,leave_plots_open=False):
#    plt.rc('text', usetex=True)

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,8)
    ax0.set_anchor('W')
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%row['track_type'],rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(1.07,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(row['sz_name'],row['track_type'], ax0, color='k', linestyle='-')
    ylim = ax0.get_ylim() #insure that all of the plots have the same zoom level in the y direction

    for j,phase_shift in enumerate(np.arange(row['phase_shift']-3*row['step'],row['phase_shift']+4*row['step'],row['step'])):
        if j>6: break #stop over plotting synthetic which can happen due to floats being used for phase shifts
        ax = fig.add_subplot(8,1,j+1, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        ax.set_ylabel(r"$\theta$=%.1f"%float(phase_shift),rotation=0,fontsize=14)
        ax.yaxis.set_label_coords(1.05,.45)
        min_proj_dis, max_proj_dis = plot_skewness_data(row,phase_shift,ax,color='k',linestyle='-',picker=leave_plots_open)
#        ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

    plot_chron_span_on_axes(row['sz_name'],fig.get_axes())

    ax0.set_xlim(max(min_proj_dis,min_syn_dis,-1000),min(max_proj_dis,max_syn_dis,1000))
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    fig.suptitle("%s - PhaseShift=%.1f - Step=%d"%(row['comp_name'],float(row['phase_shift']),int(row['step'])),fontsize=16)
    plot_scale_bars(ax)

    if leave_plots_open:
        plt.show()
    else:
        out_file = os.path.join(row['results_dir'],'deskew_plots',"%s_%s_%s.png"%(str(row['comp_name']),str(round(row['phase_shift'],3)),str(row['step'])))
        print("saving deskew file for %s to %s"%(row['comp_name'],out_file))
        fig.savefig(out_file)
#        with open(out_file,'w+') as fout:
#            fig.canvas.print_figure(fout)
        plt.close(fig)

def plot_skewnesses(deskew_path,leave_plots_open=False):
    deskew_df = pd.read_csv(deskew_path,sep='\t')

    row = deskew_df.iloc[0]

    comp_types = deskew_df['track_type'].tolist()

    if 'aero' in comp_types:
        if 'ship' in comp_types:
            ship_or_aero='both'
        else:
            ship_or_aero='aero'
    else: ship_or_aero='ship'

    check_generate_synth(row['sz_name'],row['age_min'],row['age_max'],ship_or_aero)

    check_dir(os.path.join(row['results_dir'],'deskew_plots'))

    for i,row in deskew_df.iterrows():
        if '#' in row['comp_name']: continue #check for and remove commented profiles
        print("starting process for %s"%row['comp_name'])
#        plot_deskew_file_skewnesses(row)
        plotting_process = Process(target = plot_deskew_file_skewnesses,args=[row],kwargs={'leave_plots_open':leave_plots_open})
        plotting_process.start()

#        if matlab:
#            if leave_plots_open: COMMAND = '''matlab -nodesktop -r "addpath('bin/MatlabCode/','raw_data'); show_skewness('%s',%s,%s,%s,%s,%s,'%s','%s','%s');"'''%(row['comp_name'],row['phase_shift'],row['step'],row['inter_lat'],row['inter_lon'],row['strike'],row['track_type'],row['data_dir'],row['results_dir'])
#            else: COMMAND = '''matlab -nodesktop -nodisplay -r "addpath('bin/MatlabCode/','raw_data'); show_skewness('%s',%s,%s,%s,%s,%s,'%s','%s','%s'); exit;"'''%(row['comp_name'],row['phase_shift'],row['step'],row['inter_lat'],row['inter_lon'],row['strike'],row['track_type'],row['data_dir'],row['results_dir'])
#            subprocess.check_call(COMMAND,shell=True)
#        else:
#            if leave_plots_open:
#                subprocess.check_call('''octave-cli --persist --eval "addpath('bin/MatlabCode/','raw_data'); show_skewness('%s',%s,%s,%s,%s,%s,'%s','%s','%s')" &'''%(row['comp_name'],row['phase_shift'],row['step'],row['inter_lat'],row['inter_lon'],row['strike'],row['track_type'],row['data_dir'],row['results_dir']),shell=True)
#            else:
#                subprocess.check_call('''octave-cli --quiet --no-window-system --eval "addpath('bin/MatlabCode/','raw_data'); show_skewness('%s',%s,%s,%s,%s,%s,'%s','%s','%s')" &'''%(row['comp_name'],row['phase_shift'],row['step'],row['inter_lat'],row['inter_lon'],row['strike'],row['track_type'],row['data_dir'],row['results_dir']),shell=True)

def plot_ridge_loc(deskew_row,ridge_loc_func,ax,**kwargs):
    data_file_path = os.path.join(deskew_row["data_dir"],deskew_row["comp_name"])
    data_df = pd.read_csv(data_file_path,names=["dist","dec_year","mag","lat","lon"],delim_whitespace=True)
    track_lon_lats = [[convert_to_0_360(lon),lat] for lon,lat in zip(data_df['lon'].tolist(),data_df['lat'].tolist())]

    ridge_lon,ridge_lat = ridge_loc_func(deskew_row['sz_name'],track_lon_lats)

    if ridge_lon==None or ridge_lat==None: print("ridge intercept values given are None"); return

    ridge_dis = Geodesic.WGS84.Inverse(float(deskew_row['inter_lat']),convert_to_0_360(deskew_row['inter_lon']),ridge_lat,ridge_lon)['s12']/1000

    ax.axvline(ridge_dis,**kwargs)

def plot_best_skewness_page(rows,results_dir,page_num,leave_plots_open=False,ridge_loc_func=None):
#    plt.rc('text', usetex=True)
#    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,1)
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%'ship',rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0],'ship',ax0, color='k', linestyle='-')
    ylim = ax0.get_ylim()

    for j,iterrow in enumerate(rows.iterrows()):
        i,row = iterrow
        ax = fig.add_subplot(8,1,j+2, sharex=ax0)
        ax.set_anchor('W')

        remove_axis_lines_and_ticks(ax)

        min_proj_dis, max_proj_dis = plot_skewness_data(row,float(row['phase_shift']),ax,color='k',linestyle='-',picker=True)

        if ridge_loc_func!=None:
            plot_ridge_loc(row,ridge_loc_func,ax,color='r',linestyle='-',alpha=1)

        ax.annotate(r"%s"%row['comp_name']+"\n"+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(row['inter_lat']),convert_to_0_360(row['inter_lon'])),xy=(-.15,.45),xycoords="axes fraction")
        ax.set_ylabel(r"$\theta$=%.1f"%float(row['phase_shift'])+"\n"+r"$e_a$=%.1f"%float(row['aei']),rotation=0,fontsize=14)
        ax.yaxis.set_label_coords(1.05,.45)
        ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

    ax = fig.add_subplot(8,1,8, sharex=ax0)
    ax.set_anchor('W')
    remove_axis_lines_and_ticks(ax)
    ax.set_ylabel("synthetic (%s)"%'aero',rotation=0,fontsize=14)
    ax.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0],'aero', ax, color='k', linestyle='-')
    ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

    plot_chron_span_on_axes(rows['sz_name'].iloc[0],fig.get_axes())

    ax0.set_xlim(-300,1200)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    if rows.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
        title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows['sz_name'].iloc[0],float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    else:
        title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    fig.suptitle(title,fontsize=16)
    plot_scale_bars(ax)

    if leave_plots_open:
        plt.show()
    else:
        out_file = os.path.join(results_dir,"page_%d.png"%page_num)
        fig.savefig(out_file)
        plt.close(fig)

def plot_best_skewnesses(deskew_path, leave_plots_open=False, best_skews_subdir="best_skews", ridge_loc_path=None):
#    dt_path,age_min,age_max,results_dir = create_matlab_datatable(deskew_path)

    deskew_df = filter_deskew_and_calc_aei(deskew_path)
    comp_types = deskew_df['track_type'].tolist()

    check_generate_synth(deskew_df['sz_name'].iloc[0],deskew_df['age_min'].iloc[0],deskew_df['age_max'].iloc[0],ship_or_aero='both')

    results_dir = os.path.join(deskew_df['results_dir'].iloc[0],best_skews_subdir)
    check_dir(results_dir)

    ridge_loc_func=None
    if ridge_loc_path!=None: ridge_loc_func = read_and_fit_ridge_data(ridge_loc_path)

    num_profiles_per_page = 6
    prev_i,page_num = 0,0
    for i in range(num_profiles_per_page,len(deskew_df.index),num_profiles_per_page):
        rows = deskew_df.iloc[prev_i:i]
        page_num = i/num_profiles_per_page
#        plot_best_skewness_page(rows, results_dir, page_num, leave_plots_open=leave_plots_open, ridge_loc_func=ridge_loc_func)
        plotting_process = Process(target = plot_best_skewness_page,args=[rows,results_dir,page_num],kwargs={'leave_plots_open':leave_plots_open,'ridge_loc_func':ridge_loc_func})
        plotting_process.start()
        prev_i = i
    rows = deskew_df.iloc[prev_i:len(deskew_df.index)]
    page_num += 1
    plotting_process = Process(target = plot_best_skewness_page,args=[rows,results_dir,page_num],kwargs={'leave_plots_open':leave_plots_open,'ridge_loc_func':ridge_loc_func})
    plotting_process.start()

#    if matlab:
#        if leave_plots_open:
#            subprocess.check_call('''matlab -nodesktop -r "addpath('bin/MatlabCode/','raw_data/'); plot_best_skewnesses('%s',%s,%s,'%s');"'''%(dt_path,age_min,age_max,results_dir),shell=True)
#        else:
#            subprocess.check_call('''matlab -nodesktop -nodisplay -r "addpath('bin/MatlabCode/','raw_data/'); plot_best_skewnesses('%s',%s,%s,'%s'); exit;"'''%(dt_path,age_min,age_max,results_dir),shell=True)
#    else:
#        if leave_plots_open:
#            subprocess.check_call('''octave-cli --persist --eval "addpath('bin/MatlabCode/','raw_data/'); plot_best_skewnesses('%s',%s,%s,'%s');"'''%(dt_path,age_min,age_max,results_dir),shell=True)
#        else:
#            subprocess.check_call('''octave-cli --eval "addpath('bin/MatlabCode/','raw_data/'); plot_best_skewnesses('%s',%s,%s,'%s');"'''%(dt_path,age_min,age_max,results_dir),shell=True)

def plot_skewness_by_spreading_zone(deskew_path,leave_plots_open=False,ridge_loc_path=None):
    deskew_df = pd.read_csv(deskew_path,sep='\t')

    tmp_deskew_path = '.tmp_deskew_file'
    for sz in deskew_df['sz_name'].drop_duplicates():
        sz_deskew_df = deskew_df[deskew_df['sz_name']==sz].copy()
        sz_deskew_df.sort_values(by="comp_name",inplace=True)
        sz_deskew_df.sort_values(by="inter_lat",ascending=False,inplace=True)
        sz_deskew_df.to_csv(tmp_deskew_path,sep='\t',index=False)
        plot_best_skewnesses(tmp_deskew_path,leave_plots_open=leave_plots_open,ridge_loc_path=ridge_loc_path,best_skews_subdir="skewness_by_spreading_zones/%s"%str(sz))

    try: os.remove(tmp_deskew_path); os.remove("%s.dt.csv"%os.path.basename(tmp_deskew_path).split('.')[0])
    except OSError: print("trouble removing temporary deskew and dt.csv files used for this function, check the code")

def plot_spreading_rate_picks_page(rows, spreading_rate_picks, results_dir, page_num, leave_plots_open=False):

    check_dir(results_dir)

    fig = plt.figure(figsize=(16, 9), facecolor='white')

    ax0 = fig.add_subplot(8,1,1)
    remove_axis_lines_and_ticks(ax0)
    ax0.set_ylabel("synthetic (%s)"%'ship',rotation=0,fontsize=14)
    ax0.yaxis.set_label_coords(-.1,.45)
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0],'ship',ax0, color='k', linestyle='-', plot_anom_spans=True)
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

        min_proj_dis, max_proj_dis = plot_skewness_data(row,float(row['phase_shift']),ax,color='k',linestyle='-',picker=True)

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
    min_syn_dis,max_syn_dis = plot_synthetic(rows['sz_name'].iloc[0],'aero', ax, color='k', linestyle='-', plot_anom_spans=True)
    ax.set_ylim(ylim) #insure that all of the plots have the same zoom level in the y direction

    ax0.set_xlim(-1000,500)
    fig.subplots_adjust(hspace=.0) #remove space between subplot axes
    if rows.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
        title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows['sz_name'].iloc[0],float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    else:
        title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
    fig.suptitle(title,fontsize=16)
    plot_scale_bars(ax)

    if leave_plots_open:
        plt.show()
    else:
        out_file = os.path.join(results_dir,"page_%d.png"%page_num)
        print("Saving as %s"%out_file)
        fig.savefig(out_file)
        plt.close(fig)

def plot_spreading_rate_picks(deskew_path,spreading_rate_picks_path,leave_plots_open=False):
    deskew_df = pd.read_csv(deskew_path,sep='\t')
    spreading_rate_picks = pd.read_csv(spreading_rate_picks_path,sep='\t',header=0,index_col=0)

    deskew_df['comp_name'] = deskew_df['comp_name'].apply(lambda x: str(x).rstrip('.Ed .Vd .lp'))
    picked_profiles = deskew_df[deskew_df['comp_name'].isin(spreading_rate_picks.columns)].drop_duplicates(subset='comp_name')

    check_generate_synth(deskew_df['sz_name'].iloc[0],deskew_df['age_min'].iloc[0],deskew_df['age_max'].iloc[0],ship_or_aero='both')
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

def plot_spreading_rate_results(sr_path,xmin=None,xmax=None,ymin=None,ymax=None,median_or_mean='mean',title="",leave_plots_open=False):
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
    if (xmin != None) and (xmax != None): plt.xlim( (xmin, xmax) )
    if (ymin != None) and (ymax != None): plt.ylim( (ymin, ymax) )
    
    i = 0
    for mean_median,std,age,n in zip(sr_df[median_or_mean].tolist(),sr_df['std'].tolist(),nested_ages,sr_df['n'].tolist()):
        #import pdb; pdb.set_trace()
        anom = sr_df.index.where(sr_df['std'] == std)[i]
        ax.add_patch(Rectangle((age[0], mean_median - (2.0*std)/np.sqrt(n)), age[1]-age[0],  (4.0*std)/np.sqrt(n),color='grey',alpha=.5,zorder=0,linewidth=0))
        ax.annotate('%s'%anom,xy=(sum(age)/len(age),mean_median + (2.0*std)/np.sqrt(n)),va='bottom',ha='center',fontsize=8)
        ax.annotate('n=%d'%n,xy=(sum(age)/len(age),mean_median - (2.0*std)/np.sqrt(n) - 1),va='top',ha='center',fontsize=8)
        
        i += 1

    plt.xlabel("Age (Myr)")
    plt.ylabel("Spreading Half Rate (km/Myr)")
    plt.title(title)

    if leave_plots_open:
        plt.show()
    else:
        out_file = os.path.join(os.path.dirname(sr_path),os.path.basename(sr_path).split('.')[0]+"_spreading_rate_%s.png"%(median_or_mean))
        fig.savefig(out_file)
        plt.close(fig)

def plot_isochron_picks(deskew_path,spreading_rate_picks_path,leave_plots_open=False):
    iso_df,average_lon,average_lat = get_lon_lat_from_plot_picks_and_deskew_file(deskew_path,spreading_rate_picks_path)
    if iso_df.empty: print("problem with input data, aborting"); return
    average_lon = convert_to_0_360(average_lon)

    # Read data files
    deskew_df = pd.read_csv(deskew_path,sep='\t')
    results_dir = deskew_df['results_dir'].iloc[0]

    # Initialize figure
    fig = plt.figure(figsize=(16,9), dpi=80)

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
        iso_line, = m.plot(lonlats[:,0], lonlats[:,1], latlon=True, label=anom, zorder=3, marker="x")
        if anom.strip('C r n') not in list(map(lambda x: x[0], chrons_info)): chrons_info.append([anom.strip('C r n'), iso_line.get_color()])

    plot_chron_info(chrons_info,m,coord_0_360=True,zorder=2,alpha=.4,linestyle='--',dashes=(5,1))

    legend = plt.legend()
    frame = legend.get_frame()
    frame.set_alpha(.7)

    if leave_plots_open:
        plt.show()
    else:
        out_file = os.path.join(results_dir,'iso_chron_plots',os.path.basename(spreading_rate_picks_path).split('.')[0]+'.png')
        check_dir(os.path.dirname(out_file))
        print("Saving to %s"%out_file)
        fig.savefig(out_file)
        plt.close(fig)


