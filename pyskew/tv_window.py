import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import wx.lib.buttons as buttons
import wx.lib.mixins.listctrl as listmix
import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
import cartopy.feature as cfeature
import cartopy.crs as ccrs
#import pyskew.plot_gravity as pg
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

class TVWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=200,geoid=Geodesic.WGS84,resolution="10m"):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Track Viewer %s"%parent.__version__,style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        try: self.center_lon = self.parent.dsk_row["inter_lon"]
        except AttributeError: self.center_lon = 180
        self.dpi=dpi
        self.grd_file = None
        self.geo_tiff_paths = []
        self.geoid=geoid
        self.resolution = resolution

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()
        self.configure()

        self.update()

    def init_UI(self):
        spacing = 10

        #----------------Build Directory and File Buttons-----------------

#        grd_sizer = wx.BoxSizer(wx.HORIZONTAL)
#        self.grd_path = wx.TextCtrl(self.panel, id=-1, size=(100,25), style=wx.TE_READONLY)
#        self.change_grd_btn = buttons.GenButton(self.panel, id=-1, label="Add Grid",size=(176, 29))
#        self.change_grd_btn.InitColours()
#        self.Bind(wx.EVT_BUTTON, self.on_change_grd_btn, self.change_grd_btn)
#        grd_sizer.Add(self.change_grd_btn, wx.ALIGN_LEFT)
#        grd_sizer.AddSpacer(20)
#        grd_sizer.Add(self.grd_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        #------------------------------------Make DropDown Box-----------------------------------------------------#

        latlon_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Gravity Window"), wx.VERTICAL)
        proj_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Projection"), wx.VERTICAL)
        down_sample_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Downsample Factor"), wx.HORIZONTAL)

        projs = ["Mercator","Mollweide","North Polar Stereographic","South Polar Stereographic"]
        self.proj_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(100, 25), value=projs[0], choices=projs, style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_proj,self.proj_box)

        self.max_lat_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.min_lat_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.max_lon_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.min_lon_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.down_sample_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))
        self.re_render_button = wx.Button(self.panel, id=wx.ID_ANY, label='Refresh Figure',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_re_render_button, self.re_render_button)

        #Projection sizer
        proj_sizer.Add(self.proj_box, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)

        #Lat-Lon Sizer
        lat_sizer = wx.BoxSizer(wx.HORIZONTAL)
        lat_sizer.AddMany([(self.min_lat_box, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                           (self.max_lat_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])
        lon_sizer = wx.BoxSizer(wx.HORIZONTAL)
        lon_sizer.AddMany([(self.min_lon_box, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing),
                           (self.max_lon_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing)])
        latlon_sizer.AddMany([(lat_sizer, 1, wx.ALIGN_TOP|wx.EXPAND, spacing),
                              (lon_sizer, 1, wx.ALIGN_BOTTOM|wx.EXPAND, spacing)])

        #Downsample sizer with downsample box and refresh button
        down_sample_sizer.AddMany([(self.re_render_button, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing),
                                   (self.down_sample_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing)])

        #Combine projection and downsample sizers
        proj_ds_sizer = wx.BoxSizer(wx.VERTICAL)
        proj_ds_sizer.AddMany([(proj_sizer, 1, wx.ALIGN_TOP|wx.EXPAND, spacing),
                               (down_sample_sizer, 1, wx.ALIGN_BOTTOM|wx.EXPAND, spacing)])

        #Combine all in final sizer
        all_txt_btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        all_txt_btn_sizer.AddMany([(proj_ds_sizer, 1, wx.ALIGN_LEFT|wx.EXPAND, spacing),
                                   (latlon_sizer, 1, wx.ALIGN_RIGHT|wx.EXPAND, spacing)])

        #-------------------------------------Make Figure----------------------------------------------------------#

        self.fig = Figure((2, 2), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.ax = self.fig.add_subplot(111)
        psk.remove_axis_lines_and_ticks(self.ax)
        self.toolbar.Hide()
        self.plot_setting = "Zoom"
        self.toolbar.zoom()
        self.canvas.Bind(wx.EVT_MIDDLE_DOWN,self.on_middle_click_plot)
        self.canvas.Bind(wx.EVT_MOTION,self.on_move_mouse_plot)
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_select_dleft_click)

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        outer_sizer.AddMany([#(grd_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND|wx.LEFT|wx.RIGHT,spacing),
                             (all_txt_btn_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND),
                             (self.canvas,10,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND)])

        self.panel.SetSizerAndFit(outer_sizer)

    def create_menu(self):
        """
        Generates Menu
        """

        self.menubar = wx.MenuBar()

        #-----------------
        # File Menu
        #-----------------

        menu_file = wx.Menu()

        m_tiff = menu_file.Append(-1, "&Add GeoTiff\tCtrl-O", "AddTiff")
        self.Bind(wx.EVT_MENU, self.on_add_tiff, m_tiff)

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Exit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.on_close_main, m_exit)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.SetMenuBar(self.menubar)

    def configure(self):
        try:
            dsk_row = self.parent.dsk_row
            TRACKSELECTED = isinstance(dsk_row,type(None))
        except: TRACKSELECTED = False
        if TRACKSELECTED: #Set up default values
            self.min_lat_box.SetValue("%.1f"%dsk_row["inter_lat"]-10)
            self.max_lat_box.SetValue("%.1f"%dsk_row["inter_lat"]+10)
            self.min_lon_box.SetValue("%.1f"%dsk_row["inter_lon"]-10)
            self.max_lon_box.SetValue("%.1f"%dsk_row["inter_lon"]+10)
        else:# TODO
            self.min_lat_box.SetValue("%.1f"%-20.)
            self.max_lat_box.SetValue("%.1f"%30.)
            self.min_lon_box.SetValue("%.1f"%190.)
            self.max_lon_box.SetValue("%.1f"%250.)
        self.down_sample_box.SetValue("%.0f"%50.)
        self.window = [None,None,None,None]
        self.down_sample_factor = None

    #########################Update UI Funcions#############################

    def update(self): #Populates Logger and makes plot
        self.make_map()
        self.plot_gravity()
        self.plot_tracks()
        self.plot_sites()

        self.canvas.draw()

    def on_close_main(self,event):
        self.parent.tvw_open=False
        self.Destroy()

    ###################Button and Dropdown Functions#########################

#    def on_change_grd_btn(self,event):
#        dlg = wx.FileDialog(
#            self, message="Choose Grid File",
#            defaultDir=self.parent.WD,
#            defaultFile="",
#            wildcard="Grid Files (*.grd,*.nc,*.ncf)|*.grd;*.nc;*.ncf|All Files (*.*)|*.*",
#            style=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST
#            )
#        if dlg.ShowModal() == wx.ID_OK:
#            self.grd_file = dlg.GetPath()
#            self.grd_path.SetValue(self.grd_file)
#            self.update()
#            dlg.Destroy()
#        else: dlg.Destroy()

    def on_select_proj(self,event):
        self.update()

    def on_re_render_button(self,event):
        self.plot_gravity()

    ###########################Figure Funcions###############################

    def on_middle_click_plot(self,event):
        if event.LeftIsDown() or event.ButtonDClick(): event.Skip(); return
        elif self.plot_setting == "Zoom":
            self.plot_setting = "Pan"
            self.toolbar.pan('off')
        elif self.plot_setting == "Pan":
            self.plot_setting = "Zoom"
            self.toolbar.zoom()
        event.Skip()

    def on_move_mouse_plot(self,event):
        try: dsk_row = self.parent.dsk_row
        except AttributeError: event.Skip(); return
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)
        lonlat = ccrs.PlateCarree().transform_point(*pos,self.proj)
        self.plot_tracer_on_self_and_parent(dsk_row,lonlat)
        self.parent.canvas.draw()
        self.canvas.draw()
        event.Skip()

    def on_select_dleft_click(self,event):
        try: dsk_row = self.parent.dsk_row
        except AttributeError: event.Skip(); return

        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        lonlat = ccrs.PlateCarree().transform_point(*pos,self.proj)
        proj_locs = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],[lonlat[0]],[lonlat[1]],dsk_row['strike'])
        self.parent.set_new_intercept(proj_locs.iloc[0]["dist"])
        self.parent.update(event)
        self.update()

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        self.center_lon = self.parent.dsk_row["inter_lon"]
        self.update()

    def on_add_tiff(self,event):
        dlg = wx.FileDialog(
            self, message="Choose GeoTiff File",
            defaultDir=self.parent.WD,
            wildcard="Tiff Files (*.tiff)|*.tiff|All Files (*.*)|*.*",
            style=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.geo_tiff_paths.append(dlg.GetPath())
            self.update()
            dlg.Destroy()
        else: dlg.Destroy()

    def make_map(self):
        #set basemap
        try: self.fig.delaxes(self.ax)
        except AttributeError: pass
        #TODO: ADD TRANSVERSE MERCATOR AT STRIKE AS OPTION
        if self.proj_box.GetValue() == 'North Polar Stereographic':
            self.proj = ccrs.NorthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
            pgeo.make_circular_ax(self.ax)
        elif self.proj_box.GetValue() == 'South Polar Stereographic':
            self.proj = ccrs.SouthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
            pgeo.make_circular_ax(self.ax)
        elif self.proj_box.GetValue() == 'Mercator':
            self.proj = ccrs.Mercator(central_longitude=self.center_lon)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
            path = mpath.Path(np.array([[0,0],[1,0],[1,1],[0,1],[0,0]]))
            self.ax.set_boundary(path, transform=self.ax.transAxes)
        elif self.proj_box.GetValue() == 'Mollweide':
            self.proj = ccrs.Mollweide(central_longitude=self.center_lon, globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
#            path = mpath.Path(np.array([[0.05,0.05],[0.95,0.05],[0.95,0.95],[0.05,0.95],[0.05,0.05]]))
#            self.ax.set_boundary(path, transform=self.fig.transFigure)
        else: self.parent.user_warning("Projection %s not supported"%str(self.proj_box.GetValue())); return

        self.ax.set_xticks(np.arange(0, 370, 10.), crs=ccrs.PlateCarree())
        self.ax.set_yticks(np.arange(-80, 90, 10.), crs=ccrs.PlateCarree())
        self.ax.tick_params(grid_linewidth=.5,grid_linestyle=":",color="k",labelsize=8)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        self.ax.xaxis.set_major_formatter(lon_formatter)
        self.ax.yaxis.set_major_formatter(lat_formatter)
#        self.ax.gridlines(color='grey', alpha=0.5, linestyle='--',linewidth=.5)
        land = cfeature.NaturalEarthFeature('physical', 'land', '110m',edgecolor="black",facecolor="",linewidth=2)
        self.ax.add_feature(land)

#        if self.proj_box.GetValue() == 'Mercator': self.ax.gridlines(color='black', alpha=0.4, linestyle='--', draw_labels=True)

    def plot_gravity(self):
#        try: TRACKSELECTED = self.parent.dsk_row!=None
#        except: TRACKSELECTED = False
#        if not TRACKSELECTED: self.parent.user_warning("No track selected, skipping gravity plotting"); return

        new_window = [float(self.min_lon_box.GetValue()),float(self.max_lon_box.GetValue()),float(self.min_lat_box.GetValue()),float(self.max_lat_box.GetValue())]
        new_down_sample_factor = float(self.down_sample_box.GetValue())
        if new_window==self.window and new_down_sample_factor==self.down_sample_factor: return
        else: self.window,self.down_sample_factor = new_window,new_down_sample_factor

        with wx.BusyInfo("Refreshing Gravity Grid",parent=self):
#            wx.Yield()
            all_lons,all_lats,all_grav = pg.get_sandwell(self.window,self.down_sample_factor)
#            self.ax.imshow(all_grav,cmap="viridis",alpha=.75,transform=ccrs.PlateCarree(),zorder=0,extent=self.window,origin="lower")
            fcm = self.ax.contourf(all_lons, all_lats, all_grav, 60, cmap="viridis", alpha=.75, transform=ccrs.PlateCarree(), zorder=0)
            self.ax.set_extent(self.window,crs=ccrs.PlateCarree())

    def plot_tracks(self):
        try:
            dsk_row = self.parent.dsk_row
            dsk_data = self.parent.deskew_df[self.parent.deskew_df["sz_name"]==dsk_row["sz_name"]]
        except AttributeError: return
#        projected_distances = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],mag_data['lon'].tolist(),mag_data['lat'].tolist(),dsk_row['strike'])
#        dis = max(abs(projected_distances["dist"]))*1000
#        geodict1 = Geodesic.WGS84.Direct(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']-90,dis)
#        geodict2 = Geodesic.WGS84.Direct(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']+90,dis)
        self.ax.plot([geodict1["lon2"],geodict2["lon2"]],[geodict1["lat2"],geodict2["lat2"]],transform=ccrs.Geodetic(),color="black",linewidth=1,linestyle='--')

        for j,(i,row) in enumerate(dsk_data.iterrows()):
            # Read in deskewed profile
            # This is hard-coded now. It will be updated to take a list of profiles in the future
            if not os.path.isfile(infile): self.parent.user_warning("Data file %s could not be found"%infile)
            dskd = utl.open_mag_file(os.path.join(row['data_dir'], row["comp_name"]))

            # Define the angle along which to project
            perp = row["strike"]-180
            lon = dskd["lon"]
            lat = dskd["lat"]
            mag = sk.phase_shift_data(dskd["mag"].tolist(),row["phase_shift"])

            # Find distance to project
            if row["track_type"] == 'ship':
                pcol = '#000000'
                scle = 0.2*1e3
            if row["track_type"] == 'aero':
                if 'Ed' in row["comp_name"]:
                    pcol = 'purple'
                else:
                    pcol = 'darkorchid'
                scle = 0.5*1e3

            # Project amplitude onto map
            mlats,mlons = [],[]
            for i in range(len(mag)):
                gdsc = self.geoid.Direct(lat[i], lon[i], perp, mag[i]*scle)
                mlons.append(gdsc['lon2'])
                mlats.append(gdsc['lat2'])

            # Plot map elements
            deskew_tracks.append(self.ax.plot(utl.convert_to_0_360(lon), lat, '--', linewidth=1.0, transform=ccrs.PlateCarree(), color=pcol, zorder=990))
            deskew_tracks.append(self.ax.plot(utl.convert_to_0_360(mlons), mlats, '-', linewidth=1.0, transform=ccrs.PlateCarree(), color=pcol, zorder=1000))
            deskew_fill.append(self.ax.fill_between(utl.convert_to_0_360(np.array(mlons)[mag>0]), np.array(mlats)[mag>0], lat[mag>0], transform=ccrs.PlateCarree(), alpha=0.5, color=pcol))

    def plot_sites(self):
        pass #should plot present site locations of tracks

    def plot_tracer_on_self_and_parent(self,dsk_row,lonlat):
        proj_locs = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],[lonlat[0]],[lonlat[1]],dsk_row['strike'])
        self.plot_tracer_point(dsk_row,proj_locs.iloc[0]["dist"],color="red",marker="o",s=10)
        self.parent.plot_tracer_point(proj_locs.iloc[0]["dist"],linestyle='--',color='red',alpha=.5)

    def plot_tracer_point(self,dsk_row,dis,**kwargs):
        try: self.point_on_track.remove()
        except (AttributeError,ValueError) as e: pass
        geodict = Geodesic.WGS84.Direct(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']-90,dis*1000)
        self.point_on_track = self.ax.scatter(geodict["lon2"],geodict["lat2"],transform=ccrs.Geodetic(),**kwargs)



