import wx, os, sys
import numpy as np
import pandas as pd
from scipy.stats import f
pymax_found = True
try: import pyrot.max as pymax
except ImportError: pymax_found = False
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from geographiclib.geodesic import Geodesic
import wx.lib.buttons as buttons
import wx.lib.mixins.listctrl as listmix
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.ticker as mticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key
import cartopy.feature as cfeature
import cartopy.crs as ccrs
sys.path.append("../scripts/")
import calc_strikes as cs
#import pyskew.plot_gravity as pg
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

class RTPWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=100,geoid=Geodesic.WGS84,resolution="50m",center_lon=0.,fontsize=8, verbose=False):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Pole Plot %s"%parent.__version__,style=default_style, size=(600*2,600*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        if not pymax_found: self.parent.user_warning("PyRot module PyMax not found, pole plot viewer will not be usable"); self.on_close_window(-1)
        else: self.parent.rtp_open=True
        self.center_lon = center_lon
        self.dpi = dpi
        self.geoid = geoid
        self.resolution = resolution
        self.fontsize = fontsize
        self.verbose = verbose
        self.poles_to_plot = []
        self.fixed_strike_se = None

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()
        self.configure()

        self.update()

    def init_UI(self):
        spacing = 10

        #------------------------------------Make DropDown Box-----------------------------------------------------#

        latlon_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Window Boundaries"), wx.VERTICAL)
        proj_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Projection"), wx.VERTICAL)
        refresh_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Refresh Figure"), wx.HORIZONTAL)

        projs = ["North Polar Stereographic","South Polar Stereographic","Orthographic"]
        self.proj_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(100, 25), value=projs[0], choices=projs, style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_proj,self.proj_box)

        self.max_lat_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.min_lat_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.max_lon_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
        self.min_lon_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(25,25))
#        self.down_sample_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))\

        self.re_render_button = wx.Button(self.panel, id=wx.ID_ANY, label='Refresh Figure',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_re_render_button, self.re_render_button)
        self.add_pole_button = wx.Button(self.panel, id=wx.ID_ANY, label='Add Pole',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_add_pole_button, self.add_pole_button)

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
        refresh_sizer.AddMany([(self.re_render_button, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing),
                               (self.add_pole_button, 1, wx.ALIGN_RIGHT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing)])

        #Combine projection and downsample sizers
        proj_ds_sizer = wx.BoxSizer(wx.VERTICAL)
        proj_ds_sizer.AddMany([(proj_sizer, 1, wx.ALIGN_TOP|wx.EXPAND, spacing),
                               (refresh_sizer, 1, wx.ALIGN_BOTTOM|wx.EXPAND, spacing)])

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
#        self.canvas.Bind(wx.EVT_MOTION,self.on_move_mouse_plot)
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_select_dleft_click)
        self.canvas.Bind(wx.EVT_RIGHT_DCLICK, self.on_select_dright_click)

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

        m_import_poles = menu_file.Append(-1, "&Import Poles From CSV", "ImportPoles")
        self.Bind(wx.EVT_MENU, self.on_import_poles, m_import_poles)

        menu_file.AppendSeparator()
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

        self.m_remove_pole = menu_edit.Append(-1, "&Remove Last Pole\tCtrl-Z", "RemovePole")
        self.Bind(wx.EVT_MENU, self.on_remove_pole, self.m_remove_pole)

        self.m_remove_all_poles = menu_edit.Append(-1, "&Remove All Poles\tCtrl-Z", "RemovePolesAll")
        self.Bind(wx.EVT_MENU, self.on_remove_all_poles, self.m_remove_all_poles)

        self.m_add_strike_unc = menu_edit.AppendCheckItem(-1, "&Add Strike Uncertainty\tCtrl-R", "CalcStrikes")
        self.Bind(wx.EVT_MENU, self.on_add_strike_unc, self.m_add_strike_unc)

        self.m_fix_strike_unc = menu_edit.AppendCheckItem(-1, "&Fix Strike Uncertainty\tCtrl-Shift-R", "FixStrikes")
        self.Bind(wx.EVT_MENU, self.on_fix_strike_unc, self.m_fix_strike_unc)

        self.m_solve_askw = menu_edit.AppendCheckItem(-1, "&Solve for Anomalous Skewness\tCtrl-A-R", "SolveAskw")
        self.Bind(wx.EVT_MENU, self.on_solve_askw, self.m_solve_askw)

        #-----------------
        # View Menu
        #-----------------

        menu_view = wx.Menu()

        self.m_show_lunes = menu_view.AppendCheckItem(-1, "&Show Lunes", "ShowLunes")
        self.m_show_lunes.Check()
        self.Bind(wx.EVT_MENU, self.on_show_lunes, self.m_show_lunes)

        self.m_show_pole = menu_view.AppendCheckItem(-1, "&Show Pole", "ShowPole")
#        self.m_show_pole.Check()
        self.Bind(wx.EVT_MENU, self.on_show_pole, self.m_show_pole)

        self.m_show_a95 = menu_view.AppendCheckItem(-1, "&Show A95", "ShowA95")
        self.m_show_a95.Check()
        self.Bind(wx.EVT_MENU, self.on_show_a95, self.m_show_a95)

        self.m_show_selected = menu_view.AppendCheckItem(-1, "&Show Selected", "ShowSelected")
        self.m_show_selected.Check()
        self.Bind(wx.EVT_MENU, self.on_show_selected, self.m_show_selected)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_view, "&View")
        self.SetMenuBar(self.menubar)

    def configure(self):
        self.min_lat_box.SetValue("%.1f"%60.)
        self.max_lat_box.SetValue("%.1f"%90.)
        self.min_lon_box.SetValue("%.1f"%-180.)
        self.max_lon_box.SetValue("%.1f"%180.)
        self.window = [None,None,None,None]

    #########################Update UI Funcions#############################

    def update(self): #Populates Logger and makes plot
        self.make_map() #Make Background Map

        if self.m_show_pole.IsChecked():

            if len(self.parent.deskew_df[self.parent.deskew_df["track_type"]=="ship"]) > 0:
                self.parent.save_max_file(".tmp.max",ship_only=True) #Save tmp max file to disk for debugging purposes and to more easily punch data into format using previous functions
                comment,header,ship_data = pymax.read_max_file(".tmp.max") #Read max file
                if len(ship_data["phs"])>2: #If more than 2 profiles (even-determined) invert ship
                    (plat,plon,pmag,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(ship_data, trial_pole=header[:3], out_path="synth_mag_gui.maxout", save_full_data_kernel=self.verbose, step=header[-1], max_steps=100, comment=comment)
                    s1_ship = np.sqrt(chisq/dof)*ship_data["phs"][0][1][1] #ship 1sigma
                else: s1_ship = 0
            else: ship_data,s1_ship = {"phs":[["none",[0.,0.]]]},0

            if len(self.parent.deskew_df[self.parent.deskew_df["track_type"]=="aero"]) > 0:
                self.parent.save_max_file(".tmp.max",aero_only=True) #Do same for aero only data
                comment,header,aero_data = pymax.read_max_file(".tmp.max")
                if len(aero_data["phs"])>2:
                    (plat,plon,pmag,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(aero_data, trial_pole=header[:3], out_path="synth_mag_gui.maxout", save_full_data_kernel=self.verbose, step=header[-1], max_steps=100, comment=comment)
                    s1_aero = np.sqrt(chisq/dof)*aero_data["phs"][0][1][1]
                else: s1_aero = 0
            else: aero_data,s1_aero = {"phs":[["none",[0.,0.]]]},0

            self.parent.save_max_file(".tmp.max",ship_1s=s1_ship,aero_1s=s1_aero) #now read all data and change s1 to match above
            comment,header,data = pymax.read_max_file(".tmp.max")
            if len(data["phs"])==0: raise ValueError("Empty Max file read during pole calculation")
#            for i in range(len(data["phs"])):
#                if len(ship_data["phs"]) > 0 and data["phs"][i][1][1]==ship_data["phs"][0][1][1]:
#                    data["phs"][i][1][1] = s1_ship
#                elif len(aero_data["phs"]) > 0 and data["phs"][i][1][1]==aero_data["phs"][0][1][1]:
#                    data["phs"][i][1][1] = s1_aero

            if self.m_solve_askw.IsChecked(): (plat,plon,pmag,askw,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data, trial_pole=header[:3], out_path="synth_mag_gui.maxout", save_full_data_kernel=self.verbose, step=header[-1], max_steps=100, comment=comment, solve_anom_skew=self.m_solve_askw.IsChecked())
            else: (plat,plon,pmag,maj_se,min_se,phi),chisq,dof = pymax.max_likelihood_pole(data, trial_pole=header[:3], out_path="synth_mag_gui.maxout", save_full_data_kernel=self.verbose, step=header[-1], max_steps=100, comment=comment, solve_anom_skew=self.m_solve_askw.IsChecked())
            if self.m_add_strike_unc.IsChecked(): #If strike unc is to be included calculate it!!!
                (maj_se,min_se,phi) = cs.calc_strikes_and_add_err(self.parent.deskew_path,max_file=".tmp.max",geoid=self.geoid,outfile=".tmp_dsk_cs",filter_by_quality=False,visualize=False,convergence_level=1e-5,solve_anom_skew=self.m_solve_askw.IsChecked(),fixed_strike_se=self.fixed_strike_se)
                os.remove(".tmp_dsk_cs")

            #write pole coordinates and 1sigmas to plot for user
            if phi<0: phi = phi+180
            elif phi>180: phi = phi%180
            if self.m_show_a95.IsChecked():
                f_factor = f.ppf(.95,2,dof)
                print(f_factor)
                maj_se,min_se = maj_se*np.sqrt(f_factor),min_se*np.sqrt(f_factor)
            if self.m_solve_askw.IsChecked(): self.ax.annotate(r"%.1f$^\circ$N, %.1f$^\circ$E"%(plat,plon)+"\n"+r"%.1f$^\circ$, %.1f$^\circ$, N%.1fE"%(maj_se,min_se,phi)+"\n"+"Anom. Skw. = %.1f"%askw+"\n"+r"$\chi^2_\nu$ = %.2f"%(chisq/dof)+"\n"+r"$1\sigma_{aero}$=%.1f"%(s1_aero)+"\n"+r"$1\sigma_{ship}$=%.1f"%(s1_ship),xy=(1-0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize,ha='right',va='top')
            else: self.ax.annotate(r"%.1f$^\circ$N, %.1f$^\circ$E"%(plat,plon)+"\n"+r"%.1f$^\circ$, %.1f$^\circ$, N%.1fE"%(maj_se,min_se,phi)+"\n"+r"$\chi^2_\nu$ = %.2f"%(chisq/dof)+"\n"+r"$1\sigma_{aero}$=%.1f"%(s1_aero)+"\n"+r"$1\sigma_{ship}$=%.1f"%(s1_ship),xy=(1-0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize,ha='right',va='top')
            #plot inverted pole
            self.ax = psk.plot_pole(plon,plat,phi,np.sqrt(chisq/dof)*maj_se,np.sqrt(chisq/dof)*min_se,m=self.ax, alpha=.5, zorder=10000)
        if self.m_show_lunes.IsChecked():
            #filter deskew_df to only data labeled "good" and plot lunes
            if self.m_solve_askw.IsChecked():
                srf,asf = self.parent.get_srf_asf()
                new_asf = lambda sr: asf(sr)+askw
                self.parent.deskew_df = sk.calc_aei(self.parent.deskew_df,srf,new_asf)
            else:
                self.parent.deskew_df = sk.calc_aei(self.parent.deskew_df,*self.parent.get_srf_asf())
            dsk_to_plot = self.parent.deskew_df[self.parent.deskew_df["quality"]=="g"]
            if self.m_show_selected.IsChecked():
                try: self.ax = psk.plot_lunes(dsk_to_plot,self.ax,idx_selected=self.parent.dsk_idx)
                except AttributeError: self.ax = psk.plot_lunes(dsk_to_plot,self.ax) #catch no selected data case
            else: self.ax = psk.plot_lunes(dsk_to_plot,self.ax)
#            os.remove(".tmp.max") #remove the deskew file on disk

        #plot any additional poles
        for pole_rec in self.poles_to_plot:
            print(pole_rec)
            self.ax = psk.plot_pole(*pole_rec[0],color=pole_rec[1],m=self.ax,zorder=1000)

        #set the map extent to match user input
        print([float(self.min_lon_box.GetValue()),float(self.max_lon_box.GetValue()),float(self.min_lat_box.GetValue()),float(self.max_lat_box.GetValue())])
        self.ax.set_extent([float(self.min_lon_box.GetValue()),float(self.max_lon_box.GetValue()),float(self.min_lat_box.GetValue()),float(self.max_lat_box.GetValue())], ccrs.PlateCarree())

        self.canvas.draw() #rerender

    def on_close_main(self,event):
        self.parent.rtp_open=False
        self.Destroy()

    ############################Menu Funcions################################

    def on_import_poles(self,event):
        dlg = wx.FileDialog(
            self, message="Choose CSV File",
            defaultDir=self.parent.WD,
            wildcard="Files (*.csv)|*.csv|All Files (*.*)|*.*",
            style=wx.FD_OPEN
            )
        if dlg.ShowModal() == wx.ID_OK:
            import_csv=dlg.GetPath()
            df_poles = pd.read_csv(import_csv,sep=None)
        dlg.Destroy()

        uniform_cols = list(map(lambda x: str(x).lower(), df_poles.columns))
        df_poles.columns = uniform_cols
        try:
            lat_col = next(filter(lambda x: x.startswith("lat"), uniform_cols))
            lon_col = next(filter(lambda x: x.startswith("lon"), uniform_cols))
            maj_col = next(filter(lambda x: x.startswith("maj"), uniform_cols))
            min_col = next(filter(lambda x: x.startswith("min"), uniform_cols))
            azi_col = next(filter(lambda x: x.startswith("azi"), uniform_cols))
        except:
            self.parent.user_warning("""Couldn't find a required column. There must be at least 5 columns.
                                        These 5 columns must have labels that start with lat, lon, maj, min, azi
                                        in any order and case insensitive. If more than one column fits these
                                        conditions then the first column is taken.""")
            return
        try: color_col = next(filter(lambda x: "color" in x, uniform_cols))
        except: color_col=None
        for i,row in df_poles.iterrows():
            if isinstance(color_col,type(None)):
                self.parent.user_warning("No Color for Pole (%.1f,%.1f), please specify"%(row[lat_col],row[lon_col]))
                cdlg = wx.ColourDialog(self)
                if cdlg.ShowModal() == wx.ID_OK:
                    color = tuple(np.array(cdlg.GetColourData().GetColour().Get())/255)
                else: color = "tab:blue"
            else: color = row[color_col]
            self.poles_to_plot.append([row[[lon_col,lat_col,azi_col,maj_col,min_col]].values,color])
        self.update()

    def on_save_plot(self,event):
        self.toolbar.save_figure()

    def on_remove_pole(self,event):
        self.poles_to_plot = self.poles_to_plot[:-1]
        self.update()

    def on_remove_all_poles(self,event):
        self.poles_to_plot = []
        self.update()

    def on_add_strike_unc(self,event):
        self.update()

    def on_fix_strike_unc(self,event):
        if self.m_fix_strike_unc.IsChecked():
            dlg = wx.TextEntryDialog(self, "Enter fixed strike 1 sigma:", "Fixed Strike Perturbation")
            if dlg.ShowModal() == wx.ID_OK:
                try: self.fixed_strike_se = float(dlg.GetValue())
                except ValueError:
                    self.parent.user_warning("Non-numeric value read as fixed strike 1 sigam, please provide a number")
                    dlg.Destroy()
                    return
            dlg.Destroy()
        else:
            self.parent.user_warning("Returning to calculated strike uncertainties")
            self.fixed_strike_se = None
        self.update()

    def on_solve_askw(self,event):
        self.update()

    def on_show_lunes(self,event):
        self.update()

    def on_show_pole(self,event):
        self.update()

    def on_show_a95(self,event):
        self.update()

    def on_show_selected(self,event):
        self.update()

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
        self.update()

    def on_add_pole_button(self,event):
        pdlg = PoleDialog(self) #run text entry dialog
        if pdlg.ShowModal() == wx.ID_OK:
            new_pole = [pdlg.lon,pdlg.lat,pdlg.phi,pdlg.a,pdlg.b]
        else: return
        pdlg.Destroy()
        cdata = wx.ColourData()
        cdata.SetChooseFull(True)
#        cdata.SetChooseAlpha(True)
        cdata.SetColour(wx.Colour(255, 0, 0, 128))
        cdlg = wx.ColourDialog(self,cdata)
        if cdlg.ShowModal() == wx.ID_OK:
            new_color = tuple(np.array(cdlg.GetColourData().GetColour().Get())/255)
        else: color = "#00FFFF88"
        if len(new_color)==3 or new_color[3]==1.: new_color =  (new_color[0],new_color[1],new_color[2],.5)
        elif len(new_color)<3:
            raise RuntimeError("If you're looking at this error in the terminal while running SynthMag GUI, you shouldn't be able to get here and something is significantly wrong with the color picker. Contact the dev on github.")
        cdlg.Destroy()
        self.poles_to_plot.append([new_pole,new_color]) #add new pole to list
        self.update() #update figure

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

#    def on_move_mouse_plot(self,event):
#        try: dsk_row = self.parent.dsk_row
#        except AttributeError: event.Skip(); return
#        pos=event.GetPosition()
#        width, height = self.canvas.get_width_height()
#        pos = [pos[0],height-pos[1]]
#        pos = self.ax.transData.inverted().transform(pos)
#        lonlat = ccrs.PlateCarree().transform_point(*pos,self.proj)
#        self.plot_tracer_on_self_and_parent(dsk_row,lonlat)
#        self.parent.canvas.draw()
#        self.canvas.draw()
#        event.Skip()

    def on_select_dleft_click(self,event): #TODO make rtp
        try: self.parent.dsk_row
        except AttributeError: event.Skip(); return

        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        plonlat = ccrs.PlateCarree().transform_point(*pos,self.proj)
        srf,asf = self.parent.get_srf_asf()
        reduced_skewness,rel_reduced_amplitude = sk.reduce_dsk_row_to_pole(self.parent.dsk_row,*plonlat,asf,srf)
        self.parent.phase_shift_box.SetValue("%.1f"%reduced_skewness)
        self.parent.deskew_df.at[self.parent.dsk_idx,'phase_shift'] = reduced_skewness
        self.parent.deskew_df.at[self.parent.dsk_idx,'rel_amp'] = rel_reduced_amplitude
        self.parent.dsk_row = self.parent.deskew_df.loc[self.parent.dsk_idx].iloc[0]

        self.parent.update(event)
        self.update()

    def on_select_dright_click(self,event): #TODO make rtp
        try: self.parent.deskew_df
        except AttributeError: event.Skip(); return

        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        plonlat = ccrs.PlateCarree().transform_point(*pos,self.proj)
        srf,asf = self.parent.get_srf_asf()
        for i,row in self.parent.deskew_df.iterrows():
            reduced_skewness,rel_reduced_amplitude = sk.reduce_dsk_row_to_pole(row,*plonlat,asf,srf)
            self.parent.deskew_df.at[i,'phase_shift'] = reduced_skewness
            self.parent.deskew_df.at[i,'rel_amp'] = rel_reduced_amplitude
        self.parent.deskew_df = sk.calc_aei(self.parent.deskew_df,srf,asf)
        try: self.parent.dsk_row = self.parent.deskew_df.loc[self.parent.dsk_idx].iloc[0]
        except (AttributeError,KeyError) as e: pass

        self.parent.update(event)
        self.update()

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        self.update()

    def make_map(self):
        #set basemap
        try: self.fig.delaxes(self.ax)
        except AttributeError: self.parent.user_warning("Unable to remove previous axis and refresh map, raise issue with Dev.")
        #TODO: ADD TRANSVERSE MERCATOR AT STRIKE AS OPTION
        if self.proj_box.GetValue() == 'North Polar Stereographic':
            self.proj = ccrs.NorthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
#            pgeo.make_circular_ax(self.ax)
        elif self.proj_box.GetValue() == 'South Polar Stereographic':
            self.proj = ccrs.SouthPolarStereo(central_longitude=self.center_lon,true_scale_latitude=None,globe=None)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
#            pgeo.make_circular_ax(self.ax)
        elif self.proj_box.GetValue() == 'Orthographic':
            self.proj = ccrs.Orthographic(central_longitude=self.center_lon)
            self.ax = self.fig.add_subplot(111,projection=self.proj)
        else: self.parent.user_warning("Projection %s not supported"%str(self.proj_box.GetValue())); return

#        self.ax.set_xticks(np.arange(0, 370, 10.), crs=ccrs.PlateCarree())
#        self.ax.set_yticks(np.arange(-80, 90, 10.), crs=ccrs.PlateCarree())
#        self.ax.tick_params(grid_linewidth=.5,grid_linestyle=":",color="k",labelsize=8)
#        lon_formatter = LongitudeFormatter(zero_direction_label=True)
#        lat_formatter = LatitudeFormatter()
#        self.ax.xaxis.set_major_formatter(lon_formatter)
#        self.ax.yaxis.set_major_formatter(lat_formatter)
        self.gl = self.ax.gridlines(color='black', alpha=1., linestyle=':',linewidth=.5, zorder=1)
        self.gl.xlocator = mticker.FixedLocator(np.arange(-180,180,20))
        self.gl.ylocator = mticker.FixedLocator(np.arange(-80,90,5))
        land = cfeature.NaturalEarthFeature('physical', 'land', self.resolution, edgecolor="black", facecolor="grey", linewidth=2,zorder=0)
        self.ax.add_feature(land)



class PoleDialog(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, wx.ID_ANY, "Pole Input", size= (650,200))
        self.panel = wx.Panel(self,wx.ID_ANY)
        self.parent = parent

        self.lblcoords = wx.StaticText(self.panel, label="Coords", pos=(20,20))
        self.tc_lat = wx.TextCtrl(self.panel, value="", pos=(110,20), size=(150,-1))
        self.tc_lon = wx.TextCtrl(self.panel, value="", pos=(110+500/3,20), size=(150,-1))

        self.lblell = wx.StaticText(self.panel, label="Unc Ell", pos=(20,60))
        self.tc_a = wx.TextCtrl(self.panel, value="", pos=(110,60), size=(150,-1))
        self.tc_b = wx.TextCtrl(self.panel, value="", pos=(110+500/3,60), size=(150,-1))
        self.tc_phi = wx.TextCtrl(self.panel, value="", pos=(110+(2*500)/3,60), size=(150,-1))

        self.saveButton = wx.Button(self.panel, id=wx.ID_OK, label="Ok", pos=(110,100))
        self.closeButton = wx.Button(self.panel, label="Cancel", pos=(210,100))
        self.saveButton.Bind(wx.EVT_BUTTON, self.SaveConnString)
        self.closeButton.Bind(wx.EVT_BUTTON, self.OnQuit)
        self.Bind(wx.EVT_CLOSE, self.OnQuit)

    def OnQuit(self,event):
        self.EndModal(wx.ID_CANCEL)

    def SaveConnString(self, event):
        try:
            self.lat = float(self.tc_lat.GetValue())
            self.lon = float(self.tc_lon.GetValue())
            self.a = float(self.tc_a.GetValue())
            self.b = float(self.tc_b.GetValue())
            self.phi = float(self.tc_phi.GetValue())
        except ValueError: self.parent.user_warning("At least one value was non-numeric"); return
        self.EndModal(wx.ID_OK)


