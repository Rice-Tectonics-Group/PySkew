import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import wx.lib.buttons as buttons
import wx.lib.mixins.listctrl  as  listmix
from netCDF4 import Dataset as netcdf_dataset
import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
from collections import OrderedDict
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import rasterio
import pyproj

class EAIWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Effective Inclination Viewer",style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi=dpi

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

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

#        proj_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Projection"), wx.VERTICAL)

#        projs = ["Mercator","Mollweide","North Polar Stereographic","South Polar Stereographic"]
#        self.proj_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(150, 25), value=projs[0], choices=projs, style=wx.CB_DROPDOWN|wx.TE_READONLY)
#        self.Bind(wx.EVT_COMBOBOX, self.on_select_proj,self.proj_box)

#        proj_sizer.Add(self.proj_box, 1, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.ALL, spacing)

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
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_select_dleft_click)

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        outer_sizer.AddMany([#(grd_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND|wx.LEFT|wx.RIGHT,spacing),
#                             (proj_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND),
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

#        m_tiff = menu_file.Append(-1, "&Add GeoTiff\tCtrl-O", "AddTiff")
#        self.Bind(wx.EVT_MENU, self.on_add_tiff, m_tiff)

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Exit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.on_close_main, m_exit)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.SetMenuBar(self.menubar)

    #########################Update UI Funcions#############################

    def update(self): #Populates Logger and makes plot
        self.ax.clear()
        self.plot_eai()
        self.canvas.draw()

    def on_close_main(self,event):
        self.parent.eai_open=False
        self.Destroy()

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

    def on_select_dleft_click(self,event):
        try: dsk_row = self.parent.dsk_row
        except AttributeError: event.Skip(); return

        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        try: dsk_df = self.parent.deskew_df
        except AttributeError: return
        min_dis,min_row = np.inf,None
        ylim = self.ax.get_ylim()
        for i,row in dsk_df.iterrows():
            dis = ((row["aei"]-pos[0])/ylim[0])**2 + ((row["inter_lat"]-pos[1])/ylim[1])**2
            if dis < min_dis:
                min_dis = dis
                min_row = row

        self.parent.track_box.SetValue(min_row["comp_name"])
        self.parent.on_select_track(event)
#        self.update()

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        self.update()

    def plot_eai(self):
        try: dsk_df = self.parent.deskew_df
        except AttributeError: return
        try: dsk_idx = self.parent.dsk_idx
        except AttributeError: dsk_idx = None
        for i,row in dsk_df.iterrows():
            if row["track_type"]=="aero":
                if "Ed.lp" in row["comp_name"]:
                    other_track = row["comp_name"].replace("Ed.lp","Vd.lp")
                elif "Hd.lp" in row["comp_name"]:
                    other_track = row["comp_name"].replace("Hd.lp","Vd.lp")
                elif "Vd.lp" in row["comp_name"]:
                    continue
#                    other_track = row["comp_name"].replace("Vd.lp","Ed.lp")
#                    if other_track not in dsk_df["comp_name"].tolist(): other_track = row["comp_name"].replace("Vd.lp","Hd.lp")
                else: self.parent.user_warning("Improperly named component files should have either Ed.lp, Hd.lp, or Vd.lp got: %s"%row["comp_name"]); return
                other_dsk_row = dsk_df[dsk_df["comp_name"]==other_track].iloc[0]
                aei = (row["aei"]+other_dsk_row["aei"])/2
            else: aei = row["aei"]
            if dsk_idx==i: marker = "s"
            else: marker = "o"
            self.ax.scatter(aei,row["inter_lat"],marker=marker,facecolor=(row["r"],row["g"],row["b"]),edgecolor="k",label=row["sz_name"])
        handles,labels = self.ax.get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        self.ax.legend(by_label.values(), by_label.keys(), fontsize=10, framealpha=.7)


