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
import matplotlib.patheffects as pe
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
from collections import OrderedDict
import pmagpy.ipmag as ipmag
import pandas as pd
import pyrot.max as pymax
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from pyskew.rtp_window import PoleDialog

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
        self.poles = []
        self.fzs = []
        self.geoid = Geodesic(6371.,0.)

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

        #-----------------

#        menu_file.AppendSeparator()
        submenu_save_plots = wx.Menu()

        m_save_plot = submenu_save_plots.Append(-1, "&Save Plot", "")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_save_plot,"save-plot")

        m_new_sub_plots = menu_file.Append(-1, "&Save Result", submenu_save_plots)

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Exit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.on_close_main, m_exit)

        #-----------------
        # Edit Menu
        #-----------------

        menu_edit = wx.Menu()

        self.m_add_pole = menu_edit.Append(-1, "&Add Pole", "")
        self.Bind(wx.EVT_MENU, self.on_add_pole, self.m_add_pole)

        self.m_add_fz = menu_edit.Append(-1, "&Add FZ Line", "")
        self.Bind(wx.EVT_MENU, self.on_add_fz, self.m_add_fz)

        #-----------------
        # View Menu
        #-----------------

        menu_view = wx.Menu()

        self.m_plot_legend = menu_view.AppendCheckItem(-1, "&Plot Legend\tCtrl-L", "PlotLeg")
        self.m_plot_legend.Check()
        self.Bind(wx.EVT_MENU, self.on_plot_legend, self.m_plot_legend)

        self.m_show_selected = menu_view.AppendCheckItem(-1, "&Show Selected\tCtrl-R", "PlotSel")
        self.m_show_selected.Check()
        self.Bind(wx.EVT_MENU, self.on_show_selected, self.m_show_selected)

        self.m_show_bad = menu_view.AppendCheckItem(-1, "&Show Unused Data\tCtrl-B", "PlotUD")
        self.m_show_bad.Check()
        self.Bind(wx.EVT_MENU, self.on_show_bad, self.m_show_bad)

        self.m_show_both_aero = menu_view.AppendCheckItem(-1, "&Show Both Components\tCtrl-A", "ShowBothComp")
        self.Bind(wx.EVT_MENU, self.on_show_both_components, self.m_show_both_aero)

        self.m_show_paleo_eq = menu_view.AppendCheckItem(-1, "&Show Paleo-Equator\tCtrl-E", "ShowPE")
        self.Bind(wx.EVT_MENU, self.on_show_paleo_eq, self.m_show_paleo_eq)

        self.m_change_fontsize = menu_view.Append(-1, "&Change fontsize", "")
        self.Bind(wx.EVT_MENU, self.on_change_fontsize, self.m_change_fontsize)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_view, "&View")
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
            except ValueError: self.parent.user_warning("Value entered was non-numeric canceling fontsize change.")
        dlg.Destroy()
        for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                     self.ax.get_xticklabels() + self.ax.get_yticklabels() + self.ax.get_legend().get_texts()):
            item.set_fontsize(self.fontsize)
        self.canvas.draw()

    def on_save_plot(self,event):
        self.toolbar.save_figure()

    def on_add_pole(self,event):
        pdlg = PoleDialog(self) #run text entry dialog
        if pdlg.ShowModal() == wx.ID_OK:
            new_pole = [pdlg.lon,pdlg.lat,pdlg.phi,pdlg.a,pdlg.b]
        else: return
        pdlg.Destroy()
        cdlg = wx.ColourDialog(self)
        if cdlg.ShowModal() == wx.ID_OK:
            new_color = np.array(cdlg.GetColourData().GetColour().Get())/255
        else: return
        cdlg.Destroy()
        self.poles.append([new_pole,tuple(new_color)]) #add new pole to list
        self.update() #update figure

    def on_add_fz(self,event):
        fzdlg = FZDialog(self) #run text entry dialog
        if fzdlg.ShowModal() == wx.ID_OK:
            new_fz = [fzdlg.name,fzdlg.lat,fzdlg.fontsize]
        else: return
        fzdlg.Destroy()
        cdlg = wx.ColourDialog(self)
        if cdlg.ShowModal() == wx.ID_OK:
            new_color = np.array(cdlg.GetColourData().GetColour().Get())/255
            new_fz.append(new_color)
        else: return
        cdlg.Destroy()
        self.fzs.append(new_fz) #add new pole to list
        self.update() #update figure

    def on_show_both_components(self,event):
        self.update()

    def on_plot_legend(self,event):
        self.update()

    def on_show_selected(self,event):
        self.update()

    def on_show_bad(self,event):
        self.update()

    def on_show_paleo_eq(self,event):
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
        try: dsk_df = self.parent.deskew_df.copy()
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
        try: dsk_df = self.parent.deskew_df.copy()
        except AttributeError: return
        try: dsk_idx = self.parent.dsk_idx
        except AttributeError: dsk_idx = None

        ###########################################Some shit Daniel did#####################################
#        # Add predicted Effective Remanent Inclination curve
#        dsk_df.sort_values(by=['inter_lat'],inplace=True)
#        srf_path = self.parent.spreading_rate_path
#        asf_path = self.parent.anomalous_skewness_path
#        srf,_ = sk.generate_spreading_rate_model(srf_path)
#        asf = sk.generate_anomalous_skewness_model(asf_path)
#        sk.create_max_file(dsk_df,srf,asf,outfile='temp.max')
#        _,_,max_file = pymax.read_max_file('temp.max')
#        pole,chisq,dof = pymax.max_likelihood_pole(max_file)
#        pred_eai = []

#        # Get unique spreading zone names
#        sz_names = dsk_df['sz_name'].unique()
#        # For each spreading zone, find all implied great-circle poles
#        iso_lat,iso_lon,iso_str = [],[],[]
#        for i,sz in enumerate(sz_names.tolist()):
#            sz_df = dsk_df[dsk_df['sz_name'] == sz]
#            gc_poles_lat,gc_poles_lon = [],[]
#            for j,row in sz_df.iterrows():
#                lat = row['inter_lat']
#                lon = row['inter_lon']
#                azi = row['strike'] - 90

#                gdsc = self.geoid.ArcDirect(lat,lon,azi,90)
#                gc_poles_lat.append(gdsc['lat2'])
#                gc_poles_lon.append(gdsc['lon2'])
#            sz_mean_pole = ipmag.fisher_mean(gc_poles_lon,gc_poles_lat)
#            start_loc = sz_df.iloc[0]
#            final_loc = sz_df.iloc[-1]

#            gdsc = self.geoid.Inverse(start_loc['inter_lat'],start_loc['inter_lon'],sz_mean_pole['inc'],sz_mean_pole['dec'])
#            start_azi = gdsc['azi2']-180
#            gdsc = self.geoid.Inverse(final_loc['inter_lat'],final_loc['inter_lon'],sz_mean_pole['inc'],sz_mean_pole['dec'])
#            final_azi = gdsc['azi2']-180
#            
#            sz_arc = np.linspace(start_azi, final_azi, num=20)
#            for i,azimuth in enumerate(sz_arc):
#                gdsc = self.geoid.ArcDirect(sz_mean_pole['inc'],sz_mean_pole['dec'],azimuth,90)
#                iso_lat.append(gdsc['lat2'])
#                iso_lon.append(gdsc['lon2'])
#                iso_str.append(gdsc['azi2'] - 90)

#        iso_df = pd.DataFrame(list(zip(iso_lon, iso_lat, iso_str)), columns=['Lon', 'Lat', 'Strike'])
##        print(iso_df)
#        for i,row in iso_df.iterrows():
#            plat,plon = np.deg2rad(pole[0]),np.deg2rad(pole[1])
#            slat,slon = np.deg2rad(row['Lat']),np.deg2rad(row['Lon'])
#            strike = np.deg2rad(row['Strike'])
#            dec = -np.arctan2(np.cos(plat)*np.sin(slon-plon),-np.sin(slat)*np.cos(plat)*np.cos(slon-plon)+np.cos(slat)*np.sin(plat))
#            gee = np.cos(plat)*np.cos(slat)*np.cos(plon-slon)+np.sin(plat)*np.sin(slat)
#            inc = np.arctan2(2*gee,np.sqrt(1-gee**2))
#            pred_eai.append(np.rad2deg(np.arctan2(np.tan(inc),np.sin(strike-dec))))
#        self.ax.plot(iso_df['Lat'], pred_eai)
        #####################################################Back to regularly scheduled programming

        # Get unique spreading zone names
        sorted_dsk_df = dsk_df.sort_values("inter_lat")
        sorted_dsk_df.reset_index(drop=True, inplace=True)
        sz_names = sorted_dsk_df['sz_name'].unique()

        for pole,color in self.poles:
            lats = np.linspace(round_near10(dsk_df["inter_lat"].min()-5),round_near10(dsk_df["inter_lat"].max()+5),1000)
            eis,s1_eis,lats_used,last_lat_break = [],[],[],None
            for i,sz in enumerate(sz_names.tolist()):
                sz_df = sorted_dsk_df[sorted_dsk_df['sz_name'] == sz]
                #decide which lats to use
                if len(self.fzs)>0 and (np.array(self.fzs)[:,1]>sz_df["inter_lat"].max()).any():
                    lat_break = next(x for x in sorted(np.array(self.fzs)[:,1]) if x > sz_df["inter_lat"].max()) #finds value of next largest value in list
                else:
                    try:
                        max_lat_idx = sz_df["inter_lat"].idxmax()
                        if sz_df["track_type"][max_lat_idx]=="aero": space = 2
                        else: space = 1
                        lat_break = (sz_df["inter_lat"].max()+sorted_dsk_df["inter_lat"][max_lat_idx+space])/2
                    except (KeyError,IndexError) as e: lat_break = lats[-1] #final sz
                if last_lat_break==None: sz_lats = lats[lats<=lat_break]
                else: sz_lats = lats[(last_lat_break<=lats) & (lats<=lat_break)]
                last_lat_break = lat_break
                lon = ((360+sz_df["inter_lon"])%360).mean()
                colats,s1_colats,decs,s1_decs = [],[],[],[]
                phi,a,b = pole[2:]
                for lat in sz_lats:
                    geodict = self.geoid.Inverse(pole[1], pole[0], lat, lon)
                    s1_colat = (a*b)/np.sqrt((b*np.cos(np.deg2rad(phi-geodict["azi1"])))**2 + (a*np.sin(np.deg2rad(phi-geodict["azi1"])))**2)
                    s1_dec = ((a*b)/np.sqrt((b*np.cos(np.deg2rad(phi-geodict["azi1"]-90)))**2 + (a*np.sin(np.deg2rad(phi-geodict["azi1"]-90)))**2))/np.sin(np.deg2rad(geodict["a12"]))
                    colats.append(geodict["a12"]); s1_colats.append(np.deg2rad(s1_colat)); decs.append(geodict["azi2"]), s1_decs.append(np.deg2rad(s1_dec))
                incs = np.arctan2(2*np.tan(np.deg2rad(90.-np.array(colats))),1.)
                dincs_dcolats = ((-2*((1/np.cos(np.deg2rad(90-np.array(colats))))**2))/(4*(np.tan(np.deg2rad(90-np.array(colats)))**2) + 1))
                var_incs = (np.array(s1_colats)**2)*(dincs_dcolats**2)
                eis += np.rad2deg(np.arctan2(np.tan(incs),np.sin(np.deg2rad(sz_df["strike"].mean() + 180 - np.array(decs))))).tolist()
                rproj_ang = np.deg2rad(sz_df["strike"].mean() + 180 - np.array(decs))
                dei_dinc = ((1/np.cos(incs))**2 * (1/np.sin(rproj_ang)))/(np.tan(incs)**2 * (1/np.sin(rproj_ang))**2 + 1)
                dei_ddec = (np.tan(incs) * (1/np.tan(rproj_ang)) * (1/np.sin(rproj_ang)))/(np.tan(incs)**2 * (1/np.sin(rproj_ang))**2 + 1)
                s1_eis += np.rad2deg(np.sqrt(var_incs*(dei_dinc**2) + (np.array(s1_decs)**2)*(dei_ddec**2))).tolist()
                lats_used += sz_lats.tolist()
            self.ax.plot(lats_used,eis,color=color,linewidth=1,zorder=1000)
            self.ax.plot(lats_used,eis,color="k",linewidth=1.5,zorder=999)
            self.ax.fill_between(lats_used,np.array(eis)-np.array(s1_eis),np.array(eis)+np.array(s1_eis),color=color,zorder=998,alpha=.3)
            if self.m_show_paleo_eq.IsChecked():
                idx_min = np.abs(eis).argmin()
                path_effects=[pe.Stroke(linewidth=1.5, foreground='k'), pe.Normal()]
                self.ax.axvline(lats_used[idx_min],color=color,linestyle="--",linewidth=1.0,zorder=998,path_effects=path_effects)
                self.ax.axhline(eis[idx_min],color=color,linestyle="--",linewidth=1.0,zorder=996,path_effects=path_effects)

        for fz in self.fzs:
            self.ax.axvline(fz[1],color=fz[-1],linestyle="--")
            self.ax.annotate(fz[0],xy=[fz[1],dsk_df["aei"].min()-10],color=fz[-1],fontsize=fz[2])

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

            if row["quality"]!="g":
                if not self.m_show_bad.IsChecked(): continue
                else: marker = "X"

            if self.m_show_selected.IsChecked() and (dsk_idx==i or other_idx==dsk_idx): self.ax.scatter(row["inter_lat"],aei,marker=marker,facecolor="None",edgecolor=(float(row["r"]),float(row["g"]),float(row["b"])))
            else: self.ax.scatter(row["inter_lat"],aei,marker=marker,facecolor=(float(row["r"]),float(row["g"]),float(row["b"])),edgecolor="k")

            #************************************** Experimental section **************************************
            #plat,plon = pole[0]*d2r,pole[1]*d2r
            #slat,slon = row['inter_lat']*d2r,row['inter_lon']*d2r
            #strike = row['strike']*d2r
            #dec = -np.arctan2(np.cos(plat)*np.sin(slon-plon),-np.sin(slat)*np.cos(plat)*np.cos(slon-plon)+np.cos(slat)*np.sin(plat))
            #gee = np.cos(plat)*np.cos(slat)*np.cos(plon-slon)+np.sin(plat)*np.sin(slat)
            #inc = np.arctan2(2*gee,np.sqrt(1-gee**2))
            #pred_eai.append(np.arctan2(np.tan(inc),np.sin(strike-dec))/d2r)
            #***************************************************************************************************

        #************************************** Experimental section **************************************
        #self.ax.plot(dsk_df['inter_lat'].drop_duplicates(), pred_eai)
        #***************************************************************************************************

        if self.m_plot_legend.IsChecked():
            for j,(r_col,g_col,b_col,sz_name) in dsk_df[["r","g","b","sz_name"]].drop_duplicates().iterrows():
                self.ax.scatter([],[],color=(r_col,g_col,b_col),label=sz_name,marker="s")
            self.ax.scatter([],[],color="grey",label="Ship Board Data",marker="o")
            if not self.m_show_both_aero.IsChecked(): self.ax.scatter([],[],color="grey",label="Aeromag Data",marker="s")
            else:
                self.ax.scatter([],[],color="grey",label="Vertical Aeromag Data",marker="^")
                self.ax.scatter([],[],color="grey",label="East Aeromag Data",marker=">")
            if self.m_show_selected.IsChecked(): self.ax.scatter([],[],edgecolor="grey",facecolor="None",label="Selected Data",marker="s")
            self.ax.legend(fontsize=self.fontsize, framealpha=.7)
#            handles,labels = self.ax.get_legend_handles_labels()
#            by_label = OrderedDict(zip(labels, handles))
#            self.ax.legend(by_label.values(), by_label.keys(), fontsize=self.fontsize, framealpha=.7)
        min_lat,max_lat = round_near10(dsk_df["inter_lat"].min()-5),round_near10(dsk_df["inter_lat"].max()+5)
        min_eai,max_eai = round_near10(dsk_df["aei"].min()-5),round_near10(dsk_df["aei"].max()+5)
        self.ax.set_yticks(np.arange(min_eai,max_eai+2,2),minor=True)
        self.ax.set_yticks(np.arange(min_eai,max_eai+10,10))
        self.ax.set_xticks(np.arange(min_lat,max_lat+2,2),minor=True)
        self.ax.set_xticks(np.arange(min_lat,max_lat+10,10))


        self.ax.set_xlabel("Present Latitude")
        self.ax.set_ylabel("Effective Remanent Inclination")

###########################################
##############Helper Functions#############
###########################################

def round_near10(x):
    return int(((x+5)//10)*10)

###########################################
##############Helper Dialogs#############
###########################################

class FZDialog(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, wx.ID_ANY, "FZ Input", size= (650,200))
        self.panel = wx.Panel(self,wx.ID_ANY)
        self.parent = parent

        self.lblName = wx.StaticText(self.panel, label="Name", pos=(20,20))
        self.tc_name = wx.TextCtrl(self.panel, value="", pos=(150,20), size=(300,-1))

        self.lbllat_fs = wx.StaticText(self.panel, label="Latitude/Fontsize", pos=(20,60))
        self.tc_lat = wx.TextCtrl(self.panel, value="", pos=(150,60), size=(150,-1))
        self.tc_fontsize = wx.TextCtrl(self.panel, value="12", pos=(150+500/2,60), size=(150,-1))

        self.saveButton = wx.Button(self.panel, id=wx.ID_OK, label="Ok", pos=(110,100))
        self.closeButton = wx.Button(self.panel, label="Cancel", pos=(210,100))
        self.saveButton.Bind(wx.EVT_BUTTON, self.SaveConnString)
        self.closeButton.Bind(wx.EVT_BUTTON, self.OnQuit)
        self.Bind(wx.EVT_CLOSE, self.OnQuit)

    def OnQuit(self,event):
        self.EndModal(wx.ID_CANCEL)

    def SaveConnString(self, event):
        try:
            self.name = self.tc_name.GetValue()
            self.lat = float(self.tc_lat.GetValue())
            self.fontsize = float(self.tc_fontsize.GetValue())
        except ValueError: self.parent.parent.user_warning("At least one value was non-numeric"); return
        self.EndModal(wx.ID_OK)
