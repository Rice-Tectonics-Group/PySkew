import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import wx.lib.mixins.listctrl  as  listmix
import matplotlib.path as mpath
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
import cartopy.feature as cfeature
import cartopy.crs as ccrs

class TVWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Track Viewer V0.1.1",style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        try: self.center_lon = self.parent.dsk_row["inter_lon"]
        except AttributeError: self.center_lon = 180
        self.dpi=dpi

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

        self.update()

    def init_UI(self):
        spacing = 10

        #------------------------------------Make DropDown Box-----------------------------------------------------#

        proj_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Projection"), wx.VERTICAL)

        projs = ["Mercator","Mollweide","North Polar Stereographic","South Polar Stereographic"]
        self.proj_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(150, 25), value=projs[0], choices=projs, style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_proj,self.proj_box)

        proj_sizer.Add(self.proj_box, 1, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.ALL, spacing)

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
        outer_sizer.AddMany([(proj_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND),
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

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Exit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.on_close_main, m_exit)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.SetMenuBar(self.menubar)

    #########################Update UI Funcions#############################

    def update(self): #Populates Logger and makes plot
        self.make_map()
        self.plot_track()

        self.canvas.draw()

    def on_close_main(self,event):
        self.parent.tvw_open=False
        self.Destroy()

    ###################Button and Dropdown Functions#########################

    def on_select_proj(self,event):
        self.update()

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

    def make_map(self):
        #set basemap
        try: self.fig.delaxes(self.ax)
        except AttributeError: pass
        if self.proj_box.GetValue() == 'North Polar Stereographic':
            self.proj = ccrs.NorthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
            self.ax.set_extent([-170,170,0,90],crs=ccrs.PlateCarree())
            pgeo.make_circular_ax(self.ax)
        elif self.proj_box.GetValue() == 'South Polar Stereographic':
            self.proj = ccrs.SouthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
            self.ax.set_extent([-170,170,-90,0],crs=ccrs.PlateCarree())
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
        else: self.parent.user_warning("Projection %s not supported"%str(self.proj_box.GetValue()))

        land = cfeature.NaturalEarthFeature('physical', 'land', '110m',edgecolor="black",facecolor="bisque")
        self.ax.add_feature(land)
#        if self.proj_box.GetValue() == 'Mercator': self.ax.gridlines(color='black', alpha=0.4, linestyle='--', draw_labels=True)
        self.ax.gridlines(color='black', alpha=0.4, linestyle='--',linewidth=.5)

    def plot_track(self):
        try: dsk_row = self.parent.dsk_row
        except AttributeError: return
        infile = os.path.join(dsk_row["data_dir"],dsk_row["comp_name"])
        if not os.path.isfile(infile): self.parent.user_warning("Data file %s could not be found"%infile)
        mag_data = utl.open_mag_file(infile)
        self.ax.plot(mag_data["lon"],mag_data["lat"],transform=ccrs.Geodetic(),color="black",linewidth=2)
        projected_distances = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],mag_data['lon'].tolist(),mag_data['lat'].tolist(),dsk_row['strike'])
        geodict1 = Geodesic.WGS84.Direct(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']-90,max(projected_distances["dist"])*1000)
        geodict2 = Geodesic.WGS84.Direct(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']+90,max(projected_distances["dist"])*1000)
#        geodict1 = Geodesic.WGS84.ArcDirect(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']-90,20)
#        geodict2 = Geodesic.WGS84.ArcDirect(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']+90,20)
        self.ax.plot([geodict1["lon2"],geodict2["lon2"]],[geodict1["lat2"],geodict2["lat2"]],transform=ccrs.Geodetic(),color="black",linewidth=1,linestyle='--')
#        self.fig.subplots_adjust(left=0.1, right=.9, bottom=0.1, top=.9)

    def plot_tracer_on_self_and_parent(self,dsk_row,lonlat):
        proj_locs = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],[lonlat[0]],[lonlat[1]],dsk_row['strike'])
        self.plot_tracer_point(dsk_row,proj_locs.iloc[0]["dist"],color="red",marker="o",s=10)
        self.parent.plot_tracer_point(proj_locs.iloc[0]["dist"],linestyle='--',color='red',alpha=.5)

    def plot_tracer_point(self,dsk_row,dis,**kwargs):
        try: self.point_on_track.remove()
        except (AttributeError,ValueError) as e: pass
        geodict = Geodesic.WGS84.DirectLine(dsk_row['inter_lat'],dsk_row['inter_lon'],dsk_row['strike']-90,-dis*1000).Position(dis*1000)
        self.point_on_track = self.ax.scatter(geodict["lon2"],geodict["lat2"],transform=ccrs.Geodetic(),**kwargs)



