import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import wx.lib.scrolledpanel
import wx.lib.buttons as buttons
import wx.lib.mixins.listctrl  as  listmix
import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
from collections import OrderedDict
import cartopy.feature as cfeature
import cartopy.crs as ccrs

class SkwLatWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=200,fontsize=6):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Skewness by Latitude %s"%parent.__version__,style=default_style)
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi=dpi
        self.fontsize=fontsize

#        self.panel = wx.Panel(self, wx.ID_ANY)
        self.scrolled_panel = wx.lib.scrolledpanel.ScrolledPanel(self, wx.ID_ANY)  # make the Panel

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

        self.scrolled_panel.SetAutoLayout(True)
        self.scrolled_panel.SetupScrolling()  # endable scrolling

        self.update()

    def init_UI(self):
        spacing,vpadding,hpadding = 0.,.01,.15

        #------------------------------------Make DropDown Box-----------------------------------------------------#

        sz_names_sizer = wx.StaticBoxSizer(wx.StaticBox(self.scrolled_panel, wx.ID_ANY, "Spreading Zone"), wx.HORIZONTAL)

        try:
            sz_names = self.parent.deskew_df["sz_name"].drop_duplicates().tolist()
            self.maximum_profiles = max([len(self.parent.deskew_df[self.parent.deskew_df==sz_name])for sz_name in sz_names])
        except AttributeError: sz_names,self.maximum_profiles = [""],6
        self.sz_names_box = wx.ComboBox(self.scrolled_panel, id=wx.ID_ANY, size=(300, 50), value=sz_names[0], choices=sz_names, style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_sz, self.sz_names_box)

        self.show_synth_button = wx.CheckBox(self.scrolled_panel, id=wx.ID_ANY, label="Plot Synthetic", size=(300,50))
        self.Bind(wx.EVT_CHECKBOX, self.on_show_synth_button, self.show_synth_button)

        sz_names_sizer.AddMany([(self.sz_names_box, 0, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_LEFT|wx.EXPAND|wx.ALL, spacing),
                                (self.show_synth_button, 0, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_RIGHT|wx.EXPAND|wx.ALL, spacing)])
#        self.panel.SetSizerAndFit(sz_names_sizer)

        #-------------------------------------Make Figure----------------------------------------------------------#

        canvas_sizer = wx.BoxSizer(wx.VERTICAL)

        self.fig = Figure((2., 1.*self.maximum_profiles), dpi=self.dpi)
        self.fig.subplots_adjust(top=1.-vpadding,right=1.-hpadding,left=hpadding,bottom=vpadding,wspace=.0,hspace=.0)
        self.canvas = FigCanvas(self.scrolled_panel, wx.ID_ANY, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Hide()
        self.plot_setting = "Zoom"
        self.toolbar.zoom()
        self.canvas.Bind(wx.EVT_MIDDLE_DOWN,self.on_middle_click_plot)
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_select_dleft_click)

        canvas_sizer.AddMany([(sz_names_sizer, 0, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.ALL, spacing),
                              (self.canvas, 1, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.EXPAND|wx.ALL, spacing)])

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        self.scrolled_panel.SetSizerAndFit(canvas_sizer)

#        outer_sizer = wx.BoxSizer(wx.VERTICAL)
#        outer_sizer.AddMany([(self.panel,1,wx.ALIGN_CENTER|wx.EXPAND),
#                            (self.scrolled_panel,2,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND)])

#        self.SetSizer(outer_sizer)
#        outer_sizer.Fit(self)

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
        self.fig.clear()
        self.plot_skewnesses_by_lat()
        self.canvas.draw()

    def on_close_main(self,event):
        self.parent.skw_lat_open=False
        self.Destroy()

    ##########################ComboBox Funcions##############################

    def on_select_sz(self,event):
        self.update()

    def on_show_synth_button(self,event):
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

    def on_select_dleft_click(self,event):
        pass
#        try: dsk_df = self.parent.deskew_df
#        except AttributeError: event.Skip(); return

#        pos=event.GetPosition()
#        width, height = self.canvas.get_width_height()
#        pos = [pos[0],height-pos[1]]
#        pos = self.ax.transData.inverted().transform(pos)

#        min_dis,min_row = np.inf,None
#        ylim = self.ax.get_ylim()
#        for i,row in dsk_df.iterrows():
#            dis = ((row["inter_lat"]-pos[0])/ylim[0])**2 + ((row["aei"]-pos[1])/ylim[1])**2
#            if dis < min_dis:
#                min_dis = dis
#                min_row = row

#        self.parent.track_box.SetValue(min_row["comp_name"])
#        self.parent.on_select_track(event)

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        pass

    def plot_skewnesses_by_lat(self,clip_on=True):

        try:
            sz_name = self.sz_names_box.GetValue()
            rows = self.parent.deskew_df[self.parent.deskew_df["sz_name"]==sz_name]
            rows.sort_values("inter_lat",inplace=True,ascending=False)
        except (AttributeError,KeyError) as e: print("Spreading Zone %s not found in deskew file"%str(self.sz_names_box.GetValue())); return
        try:
            xlims = self.parent.ax.get_xlim()
            ylims = self.parent.ax.get_ylim()
        except AttributeError: xlims,ylims = (-300,300),(-150,150)
        axs = self.fig.subplots(self.maximum_profiles,1,sharex=True,sharey=True)
#        for ax in axs:
#            ax.set_facecolor("grey")
        axs = axs[:len(rows)]

        for j,(ax,(i,row)) in enumerate(zip(axs,rows.iterrows())):
            print(j,row["comp_name"],xlims,ylims)
            ax.set_anchor('W')

            psk.remove_axis_lines_and_ticks(ax)

            min_proj_dis, max_proj_dis = psk.plot_skewness_data(row,float(row['phase_shift']),ax,picker=True, clip_on=clip_on, xlims=xlims, flip=True)

            ax.annotate(r"%s"%row['comp_name']+"\n"+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(row['inter_lat']),utl.convert_to_0_360(row['inter_lon'])),xy=(-.215,.5),xycoords="axes fraction",fontsize=self.fontsize,va="center",ha="left")
            ax.annotate(r"$\theta$=%.1f"%float(row['phase_shift'])+"\n"+r"$e_a$=%.1f"%float(row['aei']),xy=(1.15,.5),xycoords="axes fraction",fontsize=self.fontsize,va="center",ha="right")
#            ax.set_ylabel(r"$\theta$=%.1f"%float(row['phase_shift'])+"\n"+r"$e_a$=%.1f"%float(row['aei']),rotation=0,fontsize=self.fontsize)
            ax.yaxis.set_label_coords(1.05,.45)
            ax.patch.set_alpha(0.0)
#            ax.format_coord = format_coord

            if self.show_synth_button.GetValue():
                try:
                    ax.plot(self.parent.dis_synth,self.parent.synth,'r-',alpha=.4,zorder=1)
                except (AttributeError,IndexError): print("No synthetic found to render in skewness by latitude window")

        scale = np.sqrt(sum(np.array(xlims)**2))
        if not scale<20 or scale>3000:
            ax.set_xlim(xlims)
            ax.set_ylim(ylims)

        if self.parent.spreading_rate_path!=None:
            psk.plot_chron_span_on_axes(sz_name,self.fig.get_axes(),rows[['age_min','age_max']].iloc[0],spreading_rate_path=self.parent.spreading_rate_path)

#        self.fig.subplots_adjust(hspace=.0) #remove space between subplot axes

        self.canvas.draw()
#        if rows.groupby(level=0).agg(lambda x: len(set(x)) == 1)['sz_name'].iloc[0]:
#            title = "%s : %.1f$^\circ$N to %.1f$^\circ$N"%(rows['sz_name'].iloc[0],float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
#        else:
#            title = "%.1f$^\circ$N to %.1f$^\circ$N"%(float(rows['inter_lat'].iloc[0]),float(rows['inter_lat'].iloc[-1]))
#        fig.suptitle(title,fontsize=16)
#        psk.plot_scale_bars(ax)



