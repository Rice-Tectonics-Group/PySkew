import wx, os, sys
import numpy as np
import pandas as pd
import wx.lib.buttons as buttons
import pyskew.plot_skewness as psk
import pyskew.skewness as sk
import pyskew.utilities as utl
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

class SynthMagGUI(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,config_file=None,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, None, title="SynthMagGUI V0.0.1",style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        #Save input variables
        self.dpi=dpi
        self.WD=os.getcwd()
        self.track = None #default no data
        self.syn_buff = .2
        self.min_age,self.max_age = 0,0
        self.timescale_path,self.deskew_path = None,None

        #make the Panel
        self.panel = wx.Panel(self,-1,size=(700,450))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

        #Configure self if there is a config file
        self.configure(config_file)
        self.update(-1)

    def init_UI(self):
        """
        Builds User Interface for the interpretation Editor
        """

        side_bar_h_space = 10
        side_bar_v_space = 10
        v_bar_v_space = 10

        #----------------Build Directory and File Buttons-----------------
        wd_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Working Directory"), wx.HORIZONTAL)
        self.dir_path = wx.TextCtrl(self.panel, id=-1, size=(600,25), style=wx.TE_READONLY)
        self.change_dir_button = buttons.GenButton(self.panel, id=-1, label="change directory",size=(-1, -1))
        self.change_dir_button.SetBackgroundColour("#F8F8FF")
        self.change_dir_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_change_dir_button, self.change_dir_button)
        wd_sizer.Add(self.change_dir_button, wx.ALIGN_LEFT)
        wd_sizer.AddSpacer(40)
        wd_sizer.Add(self.dir_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        ts_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Timescale File"), wx.HORIZONTAL)
        self.ts_path = wx.TextCtrl(self.panel, id=-1, size=(600,25), style=wx.TE_READONLY)
        self.change_ts_button = buttons.GenButton(self.panel, id=-1, label="change timescale",size=(-1, -1))
        self.change_ts_button.SetBackgroundColour("#F8F8FF")
        self.change_ts_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_change_ts_button, self.change_ts_button)
        ts_sizer.Add(self.change_ts_button, wx.ALIGN_LEFT)
        ts_sizer.AddSpacer(40)
        ts_sizer.Add(self.ts_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        dsk_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Deskew File"), wx.HORIZONTAL)
        self.dsk_path = wx.TextCtrl(self.panel, id=-1, size=(600,25), style=wx.TE_READONLY)
        self.change_dsk_button = buttons.GenButton(self.panel, id=-1, label="change deskew file",size=(-1, -1))
        self.change_dsk_button.SetBackgroundColour("#F8F8FF")
        self.change_dsk_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_change_dsk_button, self.change_dsk_button)
        dsk_sizer.Add(self.change_dsk_button, wx.ALIGN_LEFT)
        dsk_sizer.AddSpacer(40)
        dsk_sizer.Add(self.dsk_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        dir_files_sizer = wx.BoxSizer(wx.VERTICAL)
        dir_files_sizer.AddMany([(wd_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space),
                          (ts_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space),
                          (dsk_sizer, 1, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space)])

        #----------------Track Selection Box-----------------

        inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        track_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Track Parameters"), wx.VERTICAL)

        self.track_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(150, 25),choices=[], style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_track,self.track_box)

        self.phase_shift_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_phase_shift, self.phase_shift_box)

        self.show_component_button = wx.CheckBox(self.panel, id=wx.ID_ANY, label="Show Other Component", size=(200,25))
        self.Bind(wx.EVT_TEXT_ENTER, self.on_show_component, self.show_component_button)

        inner_sizer.AddMany([(self.track_box, 3, wx.ALIGN_CENTER|wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.TOP,side_bar_h_space),
                             (self.phase_shift_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND|wx.TOP,side_bar_h_space)])

        track_sizer.AddMany([(inner_sizer, 1, wx.ALIGN_CENTER|wx.EXPAND),
                             (self.show_component_button, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,side_bar_h_space)])

        #----------------Age bounds sand Spreading Rate Box-----------------

        age_sr_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Age Range and Spreading Rate"), wx.HORIZONTAL)
        self.age_min_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(50, 25),choices=[], style=wx.CB_DROPDOWN|wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_age_min,self.age_min_box)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_age_min, self.age_min_box)

        self.age_max_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(50, 25),choices=[], style=wx.CB_DROPDOWN|wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_age_max,self.age_max_box)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_age_max, self.age_max_box)

        self.sr_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_sr, self.sr_box)

        age_sr_sizer.AddMany([(self.sr_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                              (self.age_min_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                              (self.age_max_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space)])

        #----------------Layer Properties and Azimuth Box-----------------

        layer_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Layer Properties and Azimuth"), wx.HORIZONTAL)
        self.layer_depth_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_layer_depth, self.layer_depth_box)

        self.layer_thickness_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_layer_thickness, self.layer_thickness_box)

        self.layer_mag_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_layer_mag, self.layer_mag_box)

        self.azi_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_azi, self.azi_box)
        layer_sizer.AddMany([(self.layer_depth_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.layer_thickness_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.layer_mag_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.azi_box, 1, wx.ALIGN_LEFT|wx.LEFT|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT,side_bar_h_space)])

        #----------------Magnetic Properties Box-----------------

        mag_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Ambient and Remanent Mag Dirs"), wx.HORIZONTAL)
        self.ad_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_ad, self.ad_box)

        self.ai_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_ai, self.ai_box)

        self.rd_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_rd, self.rd_box)

        self.ri_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_ri, self.ri_box)
        mag_sizer.AddMany([(self.ad_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.ai_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.rd_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.ri_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space)])

        #----------------Buffer Buttons Box-----------------

        zero_btn_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Fix Sythetic Endpoints"), wx.VERTICAL)
        self.zero_start_button = wx.CheckBox(self.panel, id=wx.ID_ANY, label="Zero Start Point", size=(200,25))
        self.Bind(wx.EVT_TEXT_ENTER, self.on_zero_start, self.zero_start_button)

        self.zero_end_button = wx.CheckBox(self.panel, id=wx.ID_ANY, label="Zero End Point", size=(200,25))
        self.Bind(wx.EVT_TEXT_ENTER, self.on_zero_end, self.zero_end_button)
        zero_btn_sizer.AddMany([(self.zero_start_button, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT,side_bar_h_space),
                             (self.zero_end_button, 1, wx.ALIGN_BOTTOM|wx.ALIGN_LEFT|wx.LEFT|wx.RIGHT,side_bar_h_space)])

        #----------------Transition Filter Length SynthBox-----------------

        filter_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Transition Width Filter and Synth Length"), wx.HORIZONTAL)
        self.twf_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_twf, self.twf_box)

        self.samp_n_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_samp_n, self.samp_n_box)

        self.samp_dis_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_READONLY)

        filter_sizer.AddMany([(self.twf_box, 2, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                              (self.samp_n_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space),
                              (self.samp_dis_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT,side_bar_h_space)])

        #----------------Update Button-----------------

        self.update_button = wx.Button(self.panel, id=wx.ID_ANY, label='Update',size=(200,25))
        self.Bind(wx.EVT_BUTTON, self.on_update_button, self.update_button)

        #----------------Update Button-----------------

        self.fig = Figure((5, 5), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
        self.ax = self.fig.add_subplot(111)
        psk.remove_axis_lines_and_ticks(self.ax)
        self.toolbar.Hide()
        self.plot_setting = "Zoom"
        self.toolbar.zoom()
        self.canvas.Bind(wx.EVT_MIDDLE_DOWN,self.on_middle_click_plot)
        self.canvas.Bind(wx.EVT_MOTION,self.on_move_mouse_plot)

        #------------------------------------Finish Building UI---------------------------------------------------

        #Finish Side Bar
        side_bar = wx.BoxSizer(wx.VERTICAL)
        side_bar.AddMany([(track_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (age_sr_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (layer_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (mag_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (zero_btn_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (filter_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space),
                          (self.update_button, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.TOP|wx.LEFT|wx.RIGHT|wx.EXPAND,side_bar_v_space)])

        fig_and_bar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        fig_and_bar_sizer.AddMany([(side_bar, 0, wx.ALIGN_TOP|wx.ALIGN_LEFT),
                                   (self.canvas, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.EXPAND)])

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        outer_sizer.AddMany([(dir_files_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.EXPAND),
                             (fig_and_bar_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_LEFT|wx.EXPAND|wx.TOP,v_bar_v_space)])

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

        submenu_save_plots = wx.Menu()

        m_save_deskew = submenu_save_plots.Append(-1, "&Save Deskew File", "")
        self.Bind(wx.EVT_MENU, self.on_save_deskew, m_save_deskew,"save-deskew")

        m_save_max_file = submenu_save_plots.Append(-1, "&Save Max File", "")
        self.Bind(wx.EVT_MENU, self.on_save_max_file, m_save_max_file,"save-max")

        m_save_plot = submenu_save_plots.Append(-1, "&Save Plot", "")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_save_plot,"save-plot")

        m_new_sub_plots = menu_file.Append(-1, "&Save Result", submenu_save_plots)

        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "&Exit\tCtrl-Q", "Exit")
        self.Bind(wx.EVT_MENU, self.on_close_main, m_exit)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.SetMenuBar(self.menubar)

    #########################Update UI Funcions#############################

    def configure(self,config_file):
        if config_file==None: #Set up default values
            self.dir_path.SetValue(self.WD)
            self.age_min_box.SetValue("%.1f"%float(self.min_age))
            self.age_max_box.SetValue("%.1f"%float(self.max_age))
            self.sr_box.SetValue("%.1f"%40.0)
            self.phase_shift_box.SetValue("%.1f"%0.0)
            self.layer_depth_box.SetValue("%.2f"%4.5)
            self.layer_thickness_box.SetValue("%.2f"%.5)
            self.layer_mag_box.SetValue("%.1f"%(.01*1e5))
            self.azi_box.SetValue("%.1f"%0.0)
            self.ad_box.SetValue("%.1f"%0.0)
            self.ai_box.SetValue("%.1f"%90.0)
            self.rd_box.SetValue("%.1f"%0.0)
            self.ri_box.SetValue("%.1f"%90.0)
            self.twf_box.SetValue("%.1f"%0.0)
            self.samp_n_box.SetValue("%d"%4096)
        else:# TODO
            sconf = pd.read_csv(config_path,sep='\t',header=0)
            self.dir_path.SetValue(self.WD)
            self.sr_box.SetValue("%.1f"%40.0)
            self.phase_shift_box.SetValue("%.1f"%0.0)
            self.layer_depth_box.SetValue("%.2f"%4.5)
            self.layer_thickness_box.SetValue("%.2f"%.5)
            self.layer_mag_box.SetValue("%.1f"%(.01*1e5))
            self.azi_box.SetValue("%.1f"%0.0)
            self.ad_box.SetValue("%.1f"%0.0)
            self.ai_box.SetValue("%.1f"%90.0)
            self.rd_box.SetValue("%.1f"%0.0)
            self.ri_box.SetValue("%.1f"%90.0)
            self.twf_box.SetValue("%.1f"%0.0)
            self.samp_n_box.SetValue("%d"%4096)

    def update(self,event):

        if self.timescale_path==None or self.deskew_path==None: return

        #Update Synthetic
        try:
            layer_depth = float(self.layer_depth_box.GetValue())
            layer_thickness = float(self.layer_thickness_box.GetValue())
            layer_mag = float(self.layer_mag_box.GetValue())
            azi = float(self.azi_box.GetValue())
            rd = float(self.rd_box.GetValue())
            ri = float(self.ri_box.GetValue())
            ad = float(self.ad_box.GetValue())
            ai = float(self.ai_box.GetValue())
            twf = float(self.twf_box.GetValue())
            spreading_rate = float(self.sr_box.GetValue())
            phase_shift = float(self.phase_shift_box.GetValue())
            syn_length = int(self.samp_n_box.GetValue())
        except ValueError: self.user_warning("At least one value is not numeric"); return
        fix_sta,fix_end = self.zero_start_button.GetValue(),self.zero_end_button.GetValue()
        try: synth = psk.make_synthetic(self.min_age,self.max_age,layer_depth,layer_thickness,layer_mag,azi,rd,ri,ad,ai,fix_sta,fix_end,twf,self.timescale_path,spreading_rate=spreading_rate,length=syn_length,buff=self.syn_buff)
        except: synth=[[0],[0],0]

        #Update Readouts on synthetic
        self.samp_dis_box.SetValue("%.2f"%float(synth[2]))

        #Update Data
        self.dsk_row["strike"] = (azi+90)%360
        self.dsk_row["phase_shift"] = phase_shift
        self.deskew_df[self.dsk_idx]["strike"] = (azi+90)%360
        self.deskew_df[self.dsk_idx]["phase_shift"] = phase_shift

        #Center Synthetic
        anom_width = abs(self.dsk_row["age_max"]*spreading_rate-self.dsk_row["age_min"]*spreading_rate)/2
        center_dis=((self.dsk_row["age_max"]-self.min_age)*spreading_rate+(self.dsk_row["age_max"]-self.min_age)*spreading_rate)/2 - anom_width
        neg_anom=((-self.dsk_row["age_min"]-self.dsk_row["age_max"]))*spreading_rate
        dis_synth = np.array(synth[1])-center_dis

        ylim,xlim = self.ax.get_ylim(),self.ax.get_xlim()
        self.ax.clear()
#        psk.remove_axis_lines_and_ticks(self.ax)
        if self.show_component_button.GetValue():
            if self.dsk_row["track_type"]=="aero":
                if "Ed.lp" in self.track:
                    other_track = self.track.replace("Ed.lp","Vd.lp")
                    other_phase = phase_shift-90
                elif "Vd.lp" in self.track:
                    other_track = self.track.replace("Vd.lp","Ed.lp")
                    other_phase = phase_shift+90
                else: self.user_warning("Improperly named component files should have either Ed.lp or Vd.lp got: %s"%self.track); return
                other_dsk_row = self.deskew_df[self.deskew_df["comp_name"]==other_track].iloc[0]
                other_dsk_row["strike"] = (azi+90)%360
                psk.plot_skewness_data(other_dsk_row,other_phase,self.ax,color='darkgreen',zorder=2,picker=True,alpha=.7)
            else: self.user_warning("Cannot show other componenet for track type: %s"%str(self.deskew_df["track_type"]))
        psk.plot_skewness_data(self.dsk_row,self.dsk_row["phase_shift"],self.ax,color='k',zorder=3,picker=True)
        self.ax.plot(dis_synth,synth[0],'r-',alpha=.4,zorder=1)
        self.ax.plot(dis_synth,np.zeros(len(dis_synth)),'k--')
        self.ax.axvspan(-anom_width,anom_width, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)
        if self.min_age<0: self.ax.axvspan(neg_anom-anom_width,neg_anom+anom_width, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)
#        psk.plot_scale_bars(self.ax,offset_of_bars = .05)
        self.ax.annotate("%s\n%s\n"%(self.dsk_row["sz_name"],self.track)+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(self.dsk_row['inter_lat']),utl.convert_to_0_360(self.dsk_row['inter_lon'])),xy=(0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5))
        if not xlim==(0.0,1.0):
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)

        self.canvas.draw()

    def update_age_boxes(self):
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0)

        self.age_min_box.Clear()
        self.age_min_box.SetItems(tdf["chron"].tolist())

        self.age_max_box.Clear()
        self.age_max_box.SetItems(tdf["chron"].tolist())

    def update_track_box(self):
        self.deskew_df = utl.open_deskew_file(self.deskew_path)
        abs_data_paths = [self.deskew_df["data_dir"][i] if self.deskew_df["data_dir"][i]==os.path.abspath(self.deskew_df["data_dir"][i]) else os.path.abspath(os.path.join(self.WD,self.deskew_df["data_dir"][i])) for i in self.deskew_df.index]
#        if any(not os.path.isfile(abs_data_paths)): abs_data_paths = [self.deskew_df["data_dir"][i] if self.deskew_df["data_dir"][i]==os.path.abspath(self.deskew_df["data_dir"][i]) else os.path.abspath(os.path.join(os.path.dirname(self.deskew_path),self.deskew_df["data_dir"][i])) for i in self.deskew_df.index]
        self.deskew_df["data_dir"] = abs_data_paths
        self.track_box.Clear()
        self.track_box.SetItems(self.deskew_df["comp_name"].tolist())

    def on_close_main(self,event):
        self.Destroy()

    ##########################Directory/File Buttons########################

    def on_change_dir_button(self,event):
        currentDirectory = os.getcwd()
        dlg = wx.DirDialog(self.panel, "Choose Your New Working Directory:", defaultPath=currentDirectory, style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON | wx.DD_CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            self.WD = dlg.GetPath()
            os.chdir(self.WD)
            self.dir_path.SetValue(self.WD)
            dlg.Destroy()
        else: dlg.Destroy()

    def on_change_ts_button(self,event):
        dlg = wx.FileDialog(
            self, message="Choose Timescale File",
            defaultDir=self.WD,
            defaultFile="timescale.txt",
            wildcard="Files (*.txt)|*.txt|All Files (*.*)|*.*",
            style=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.timescale_path = dlg.GetPath()
            self.ts_path.SetValue(self.timescale_path)
            self.update_age_boxes()
            dlg.Destroy()
        else: dlg.Destroy()

    def on_change_dsk_button(self,event):
        dlg = wx.FileDialog(
            self, message="Choose Deskew File",
            defaultDir=self.WD,
            defaultFile="chron*.txt",
            wildcard="Files (*.deskew)|*.deskew|All Files (*.*)|*.*",
            style=wx.FD_OPEN|wx.FD_FILE_MUST_EXIST
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.deskew_path = dlg.GetPath()
            self.dsk_path.SetValue(self.deskew_path)
            self.update_track_box()
            dlg.Destroy()
        else: dlg.Destroy()

    ##########################Drop Down Boxes################################

    def on_select_track(self,event):
        self.track = self.track_box.GetValue()
        self.dsk_row = self.deskew_df[self.deskew_df["comp_name"]==self.track].iloc[0]
        self.dsk_idx = self.deskew_df.index[self.deskew_df["comp_name"]==self.track]
        try:
            if "strike" in self.dsk_row and not np.isnan(float(self.dsk_row["strike"])):
                self.azi_box.SetValue("%.1f"%(float(self.dsk_row["strike"])-90))
            if "phase_shift" in self.dsk_row and not np.isnan(float(self.dsk_row["phase_shift"])):
                self.phase_shift_box.SetValue("%.1f"%float(self.dsk_row["phase_shift"]))
        except TypeError: self.user_warning("Invalid Strike or Phase Shift in deskew file for %s"%self.track)

    def on_select_age_min(self,event):
        min_chron = self.age_min_box.GetValue()
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        self.min_age = tdf.loc[min_chron]["top"]
        if self.min_age>self.max_age: self.min_age,self.max_age = self.swap(self.min_age,self.max_age)

    def on_select_age_max(self,event):
        max_chron = self.age_max_box.GetValue()
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        self.max_age = tdf.loc[max_chron]["base"]
        if self.min_age>self.max_age: self.min_age,self.max_age = self.swap(self.min_age,self.max_age)

    ##########################Text Entry Boxes################################

    def on_enter_age_min(self,event):
        min_value = self.age_min_box.GetValue()
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        if min_value in tdf.index:
            self.min_age = tdf.loc[min_value]["top"]
        else:
            try: self.min_age = float(min_value)
            except ValueError: self.user_warning("%s is not a valid minimum age value"%str(min_value))
        if self.min_age>self.max_age: self.min_age,self.max_age = self.swap(self.min_age,self.max_age)
        self.update(event)

    def on_enter_age_max(self,event):
        max_value = self.age_max_box.GetValue()
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        if max_value in tdf.index:
            self.max_age = tdf.loc[max_value]["base"]
        else:
            try: self.max_age = float(max_value)
            except ValueError: self.user_warning("%s is not a valid minimum age value"%str(max_value))
        if self.min_age>self.max_age: self.min_age,self.max_age = self.swap(self.min_age,self.max_age)
        self.update(event)

    def on_enter_sr(self,event):
        self.update(event)

    def on_enter_phase_shift(self,event):
        self.update(event)

    def on_enter_layer_depth(self,event):
        self.update(event)

    def on_enter_layer_thickness(self,event):
        self.update(event)

    def on_enter_layer_mag(self,event):
        self.update(event)

    def on_enter_azi(self,event):
        self.update(event)

    def on_enter_ad(self,event):
        self.update(event)

    def on_enter_ai(self,event):
        self.update(event)

    def on_enter_rd(self,event):
        self.update(event)

    def on_enter_ri(self,event):
        self.update(event)

    def on_enter_twf(self,event):
        self.update(event)

    def on_enter_samp_n(self,event):
        self.update(event)

    ##########################Radio Buttons################################

    def on_zero_start(self,event):
        pass

    def on_zero_end(self,event):
        pass

    def on_show_component(self,event):
        pass

    ##########################Buttons!!!!!!################################

    def on_update_button(self,event):
        self.update(event)
#        utl.run_in_parallel(self.update,args=[event])

    ##########################Plot Functions################################

    def on_middle_click_plot(self,event):
        if event.LeftIsDown() or event.ButtonDClick():
            return
        elif self.plot_setting == "Zoom":
            self.plot_setting = "Pan"
            self.toolbar.pan('off')
        elif self.plot_setting == "Pan":
            self.plot_setting = "Zoom"
            self.toolbar.zoom()

    def on_move_mouse_plot(self,event):
        pos=event.GetPosition()
        try: self.point_annotation.remove()
        except (AttributeError,ValueError) as e: pass
        pos = self.ax.transData.inverted().transform(pos)
        self.point_annotation = self.ax.annotate("x = %.2f\ny = %.2f"%(float(pos[0]),float(pos[1])),xy=(1-0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5))
        self.canvas.draw()
        event.Skip()

    ##########################Plot Functions################################

    def on_save_deskew(self,event):
        dlg = wx.FileDialog(
            self, message="Save Deskew File",
            defaultDir=self.WD,
            defaultFile=os.path.basename(self.deskew_path),
            wildcard="Files (*.deskew)|*.deskew|All Files (*.*)|*.*",
            style=wx.FD_SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            outfile = dlg.GetPath()
            self.deskew_df.to_csv(outfile,sep="\t",index=False)
        dlg.Destroy()

    def on_save_plot(self,event):
        try: self.point_annotation.remove()
        except (AttributeError,ValueError) as e: pass
        self.canvas.draw()
        self.toolbar.save_figure()

    def on_save_max_file(self,event):
        avg_age = (self.deskew_df.iloc[0]["age_min"]+self.deskew_df.iloc[0]["age_max"])/2
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0)
        chron = tdf[(tdf["top"]<=avg_age) & (tdf["base"]>=avg_age)].iloc[0]["chron"]

        dlg = wx.FileDialog(
            self, message="Save Max File",
            defaultDir=self.WD,
            defaultFile="%s.max"%chron,
            wildcard="Files (*.max)|*.max|All Files (*.*)|*.*",
            style=wx.FD_SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            outfile = dlg.GetPath()
            self.deskew_df.to_csv(".tmp.deskew",sep="\t",index=False)
            sk.create_maxtab_file(".tmp.deskew",chron,outfile=outfile)
        dlg.Destroy()


    ##########################Utility Dialogs and Functions################

    def user_warning(self, message, caption = 'Warning!'):
        """
        Shows a dialog that warns the user about some action

        Parameters
        ----------
        message : message to display to user
        caption : title for dialog (default: "Warning!")

        Returns
        -------
        continue_bool : True or False
        """
        dlg = wx.MessageDialog(self, message, caption, wx.OK | wx.CANCEL | wx.ICON_WARNING)
        if dlg.ShowModal() == wx.ID_OK:
            continue_bool = True
        else:
            continue_bool = False
        dlg.Destroy()
        return continue_bool

    def swap(self,a,b):
        return b,a




def main():
    #Parse system arguments
    if "-h" in sys.argv:
        help(__name__)
        sys.exit()
    config_file = None
    if '-f' in sys.argv:
        f_index = sys.argv.index('-f')
        if len(sys.argv)>=f_index+2:
            config_file = sys.argv[f_index+1]
    dpi = 200
    if '-dpi' in sys.argv:
        dpi_index = sys.argv.index('-dpi')
        if len(sys.argv)>=f_index+2:
            dpi = sys.argv[dpi_index+1]

    #Start GUI
    app = wx.App()
    app.frame = SynthMagGUI(config_file=config_file,dpi=dpi)
    app.frame.Center()
    app.frame.Show()
    app.MainLoop()

if __name__=="__main__":
    main()
