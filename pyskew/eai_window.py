import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
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

class EAIWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=200,fontsize=8):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Effective Inclination Viewer %s"%parent.__version__,style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi=dpi
        self.fontsize=fontsize

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

        self.m_plot_legend = menu_file.AppendCheckItem(-1, "&Plot Legend\tCtrl-L", "PlotLeg")
        self.m_plot_legend.Check()
        self.Bind(wx.EVT_MENU, self.on_plot_legend, self.m_plot_legend)

        self.m_show_both_aero = menu_file.AppendCheckItem(-1, "&Show Both Components\tCtrl-A", "ShowBothComp")
        self.Bind(wx.EVT_MENU, self.on_show_both_components, self.m_show_both_aero)

        self.m_change_fontsize = menu_file.Append(-1, "&Change fontsize", "")
        self.Bind(wx.EVT_MENU, self.on_change_fontsize, self.m_change_fontsize)

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

    ############################Menu Funcions################################

    def on_change_fontsize(self,event):
        dlg = wx.TextEntryDialog(self, "Enter Fontsize", caption="Edit Fontsize",
                    value=str(self.fontsize), style=wx.TextEntryDialogStyle)
        if dlg.ShowModal() == wx.ID_OK:
            try:
                self.fontsize = int(dlg.GetValue())
                mpl.rcParams.update({'font.size': self.fontsize})
            except ValueError: self.user_warning("Value entered was non-numeric canceling fontsize change.")
        dlg.Destroy()
        for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                     self.ax.get_xticklabels() + self.ax.get_yticklabels() + self.ax.get_legend().get_texts()):
            item.set_fontsize(self.fontsize)
        self.canvas.draw()

    def on_show_both_components(self,event):
        self.update()

    def on_plot_legend(self,event):
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
        try: dsk_df = self.parent.deskew_df
        except AttributeError: event.Skip(); return

        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        min_dis,min_row = np.inf,None
        ylim = self.ax.get_ylim()
        for i,row in dsk_df.iterrows():
            dis = ((row["inter_lat"]-pos[0])/ylim[0])**2 + ((row["aei"]-pos[1])/ylim[1])**2
            if dis < min_dis:
                min_dis = dis
                min_row = row

        self.parent.track_box.SetValue(min_row["comp_name"])
        self.parent.on_select_track(event)

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        pass

    def plot_eai(self):
        try: dsk_df = self.parent.deskew_df
        except AttributeError: return
        try: dsk_idx = self.parent.dsk_idx
        except AttributeError: dsk_idx = None

        for i,row in dsk_df.iterrows():

            other_idx = np.nan
            if row["track_type"]=="aero":
                if "Ed.lp" in row["comp_name"]:
                    if self.m_show_both_aero.IsChecked(): marker = ">"
                    else: marker = "s"
                    other_track = row["comp_name"].replace("Ed.lp","Vd.lp")
                elif "Hd.lp" in row["comp_name"]:
                    if self.m_show_both_aero.IsChecked(): marker = ">"
                    else: marker = "s"
                    other_track = row["comp_name"].replace("Hd.lp","Vd.lp")
                elif "Vd.lp" in row["comp_name"]:
                    if self.m_show_both_aero.IsChecked(): marker = "^"
                    else: continue
#                    other_track = row["comp_name"].replace("Vd.lp","Ed.lp")
#                    if other_track not in dsk_df["comp_name"].tolist(): other_track = row["comp_name"].replace("Vd.lp","Hd.lp")
                else: self.parent.user_warning("Improperly named component files should have either Ed.lp, Hd.lp, or Vd.lp got: %s"%row["comp_name"]); return
                if not self.m_show_both_aero.IsChecked():
                    other_dsk_row = dsk_df[dsk_df["comp_name"]==other_track].iloc[0]
                    other_idx = other_dsk_row.name
                    aei = (row["aei"]+other_dsk_row["aei"])/2
                else: aei = row["aei"]
            else: aei = row["aei"]; marker = "o"

            if row["quality"]!="g": marker = "X"

            if dsk_idx==i or other_idx==dsk_idx: self.ax.scatter(row["inter_lat"],aei,marker=marker,facecolor="None",edgecolor=(float(row["r"]),float(row["g"]),float(row["b"])))
            else: self.ax.scatter(row["inter_lat"],aei,marker=marker,facecolor=(float(row["r"]),float(row["g"]),float(row["b"])),edgecolor="k")

        if self.m_plot_legend.IsChecked():
            for j,(r_col,g_col,b_col,sz_name) in dsk_df[["r","g","b","sz_name"]].drop_duplicates().iterrows():
                self.ax.scatter([],[],color=(r_col,g_col,b_col),label=sz_name,marker="s")
            self.ax.scatter([],[],color="grey",label="Ship Board Data",marker="o")
            if not self.m_show_both_aero.IsChecked(): self.ax.scatter([],[],color="grey",label="Aeromag Data",marker="s")
            else:
                self.ax.scatter([],[],color="grey",label="Vertical Aeromag Data",marker="^")
                self.ax.scatter([],[],color="grey",label="East Aeromag Data",marker=">")
            self.ax.scatter([],[],edgecolor="grey",facecolor="None",label="Selected Data",marker="s")
            self.ax.legend(fontsize=self.fontsize, framealpha=.7)
#            handles,labels = self.ax.get_legend_handles_labels()
#            by_label = OrderedDict(zip(labels, handles))
#            self.ax.legend(by_label.values(), by_label.keys(), fontsize=self.fontsize, framealpha=.7)

        self.ax.set_xlabel("Present Latitude")
        self.ax.set_ylabel("Effective Remanent Inclination")


