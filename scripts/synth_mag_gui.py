import wx, os, sys, pdb
import numpy as np
import pandas as pd
import wx.lib.buttons as buttons
import pyskew.plot_skewness as psk
import pyskew.skewness as sk
from pyskew.srm_window import SRMWindow
from pyskew.tv_window import TVWindow
import pyskew.utilities as utl
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

class SynthMagGUI(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,config_file=None,fontsize=8,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, None, title="SynthMagGUI V0.1.1",style=default_style, size=(400*3,300*3))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        #Save input variables
        self.dpi=dpi
        self.WD=os.getcwd()
        self.track = None #default no data
        self.syn_buff = .2
        self.min_age,self.max_age = 0,0
        self.spreading_rate_path,self.anomalous_skewness_path,self.timescale_path,self.deskew_path = None,None,None,None
        self.srmw_open,self.tvw_open = False,False

        #make the Panel
        self.panel = wx.Panel(self,-1,size=(1200,900))
        self.fontsize = fontsize

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
        # (176, 29) is the size of the largest button this way everything is nicely spaced
        wd_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.dir_path = wx.TextCtrl(self.panel, id=-1, size=(600,25), style=wx.TE_READONLY)
        self.change_dir_button = buttons.GenButton(self.panel, id=-1, label="change directory",size=(176, 29))
#        self.change_dir_button.SetBackgroundColour("#F8F8FF")
        self.change_dir_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_change_dir_button, self.change_dir_button)
        wd_sizer.Add(self.change_dir_button, wx.ALIGN_LEFT)
        wd_sizer.AddSpacer(20)
        wd_sizer.Add(self.dir_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        ts_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.ts_path = wx.TextCtrl(self.panel, id=-1, size=(300,25), style=wx.TE_READONLY)
        self.change_ts_button = buttons.GenButton(self.panel, id=-1, label="change timescale",size=(176, 29))
#        self.change_ts_button.SetBackgroundColour("#F8F8FF")
        self.change_ts_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_open_ts_button, self.change_ts_button)
        ts_sizer.Add(self.change_ts_button, wx.ALIGN_LEFT)
        ts_sizer.AddSpacer(20)
        ts_sizer.Add(self.ts_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        dsk_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.dsk_path = wx.TextCtrl(self.panel, id=-1, size=(300,25), style=wx.TE_READONLY)
        self.change_dsk_button = buttons.GenButton(self.panel, id=-1, label="change deskew file",size=(176, 29))
#        self.change_dsk_button.SetBackgroundColour("#F8F8FF")
        self.change_dsk_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_open_dsk_button, self.change_dsk_button)
        dsk_sizer.Add(self.change_dsk_button, wx.ALIGN_LEFT)
        dsk_sizer.AddSpacer(20)
        dsk_sizer.Add(self.dsk_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        sr_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sr_path = wx.TextCtrl(self.panel, id=-1, size=(300,25), style=wx.TE_READONLY)
        self.change_sr_button = buttons.GenButton(self.panel, id=-1, label="change spreading rate",size=(176, 29))
#        self.change_sr_button.SetBackgroundColour("#F8F8FF")
        self.change_sr_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_open_sr_file, self.change_sr_button)
        sr_sizer.Add(self.change_sr_button, wx.ALIGN_LEFT)
        sr_sizer.AddSpacer(20)
        sr_sizer.Add(self.sr_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        as_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.as_path = wx.TextCtrl(self.panel, id=-1, size=(300,25), style=wx.TE_READONLY)
        self.change_as_button = buttons.GenButton(self.panel, id=-1, label="change anomalous skewness",size=(176, 29))
#        self.change_as_button.SetBackgroundColour("#F8F8FF")
        self.change_as_button.InitColours()
        self.Bind(wx.EVT_BUTTON, self.on_open_as_file, self.change_as_button)
        as_sizer.Add(self.change_as_button, wx.ALIGN_LEFT)
        as_sizer.AddSpacer(20)
        as_sizer.Add(self.as_path,wx.ALIGN_CENTER_VERTICAL|wx.EXPAND)

        indir1_files_sizer = wx.BoxSizer(wx.HORIZONTAL)
        indir1_files_sizer.AddMany([(ts_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.RIGHT, v_bar_v_space),
                                    (sr_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND)])

        indir2_files_sizer = wx.BoxSizer(wx.HORIZONTAL)
        indir2_files_sizer.AddMany([(dsk_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.RIGHT, v_bar_v_space),
                                    (as_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND)])

        dir_files_sizer = wx.BoxSizer(wx.VERTICAL)
        dir_files_sizer.AddMany([(wd_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space),
                          (indir1_files_sizer, 1, wx.ALIGN_TOP|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space),
                          (indir2_files_sizer, 1, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER|wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT, v_bar_v_space)])

        #----------------Track Selection Box-----------------

        inner_sizer = wx.BoxSizer(wx.HORIZONTAL)
        track_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Track Parameters"), wx.VERTICAL)

        self.track_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(150, 25),choices=[], style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_track,self.track_box)

        self.phase_shift_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_phase_shift, self.phase_shift_box)

        self.show_component_button = wx.CheckBox(self.panel, id=wx.ID_ANY, label="Show Other Component", size=(200,25))
        self.Bind(wx.EVT_TEXT_ENTER, self.on_show_component, self.show_component_button)

#        self.show_adjacent_button = wx.CheckBox(self.panel, id=wx.ID_ANY, label="Show Adjacent Components", size=(200,25))
#        self.Bind(wx.EVT_TEXT_ENTER, self.on_show_adjacent, self.show_adjacent_button)

        inner_sizer.AddMany([(self.track_box, 3, wx.ALIGN_CENTER|wx.EXPAND|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.TOP,side_bar_h_space),
                             (self.phase_shift_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND|wx.TOP,side_bar_h_space)])

        track_sizer.AddMany([(inner_sizer, 1, wx.ALIGN_CENTER|wx.EXPAND),
                             (self.show_component_button, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,side_bar_h_space)])
#                             (self.show_adjacent_button, 1, wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,side_bar_h_space)])

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
        self.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_select_dleft_click)

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

        m_open_ts_file = menu_file.Append(-1, "&Open Timescale File\tCtrl-T", "open-ts")
        self.Bind(wx.EVT_MENU, self.on_open_ts_button, m_open_ts_file)

        m_open_dskew_file = menu_file.Append(-1, "&Open Deskew File\tCtrl-D", "open-deskew")
        self.Bind(wx.EVT_MENU, self.on_open_dsk_button, m_open_dskew_file)

        m_open_sr_file = menu_file.Append(-1, "&Open Spreading Rate File\tCtrl-R", "open-sr")
        self.Bind(wx.EVT_MENU, self.on_open_sr_file, m_open_sr_file)

        m_open_as_file = menu_file.Append(-1, "&Open Anomalous Skewness File\tCtrl-M", "open-as")
        self.Bind(wx.EVT_MENU, self.on_open_as_file, m_open_as_file)

        menu_file.AppendSeparator()
        submenu_save_plots = wx.Menu()

        m_save_deskew = submenu_save_plots.Append(-1, "&Save Deskew File", "")
        self.Bind(wx.EVT_MENU, self.on_save_deskew, m_save_deskew,"save-deskew")

        m_save_max_file = submenu_save_plots.Append(-1, "&Save Max File", "")
        self.Bind(wx.EVT_MENU, self.on_save_max_file, m_save_max_file,"save-max")

        m_save_maxtab_file = submenu_save_plots.Append(-1, "&Save Maxtab File", "")
        self.Bind(wx.EVT_MENU, self.on_save_maxtab_file, m_save_maxtab_file,"save-maxtab")

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

        self.m_use_sr_model = menu_edit.AppendCheckItem(-1, "&Toggle Spreading Rate Model\tCtrl-Shift-R", "")
        self.Bind(wx.EVT_MENU, self.on_use_sr_model, self.m_use_sr_model)

        self.m_use_as_model = menu_edit.AppendCheckItem(-1, "&Toggle Anomalous Skewness Correction\tCtrl-Shift-M", "")
        self.Bind(wx.EVT_MENU, self.on_use_as_model, self.m_use_as_model)

        #-----------------
        # View Menu
        #-----------------

        menu_view = wx.Menu()

        self.m_show_anoms = menu_view.AppendCheckItem(-1, "&Show Anomaly Names\tCtrl-A", "")
        self.Bind(wx.EVT_MENU, self.on_show_anoms, self.m_show_anoms)

        self.m_show_major_anoms = menu_view.AppendCheckItem(-1, "&Show Major Anomaly Names Only\tCtrl-Shift-A", "")
        self.Bind(wx.EVT_MENU, self.on_show_major_anoms, self.m_show_major_anoms)

        #-----------------
        # Tools Menu
        #-----------------

        menu_tools = wx.Menu()

        self.m_srm_edit = menu_tools.Append(-1, "&Spreading Rate Model\tAlt-R", "")
        self.Bind(wx.EVT_MENU, self.on_srm_edit, self.m_srm_edit)

        self.m_tv = menu_tools.Append(-1, "&Track Viewer\tAlt-V", "")
        self.Bind(wx.EVT_MENU, self.on_tv, self.m_tv)

        #-----------------
        # Help Menu
        #-----------------

        menu_help = wx.Menu()

        m_open_debug = menu_help.AppendCheckItem(-1, "&Open Debugger\tCtrl-Shift-D", "")
        self.Bind(wx.EVT_MENU, self.on_open_debug, m_open_debug)

        #-----------------

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_view, "&View")
        self.menubar.Append(menu_tools, "&Tools")
        self.menubar.Append(menu_help, "&Help")
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

        if self.timescale_path==None: return

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
        try:
            if self.m_use_sr_model.IsChecked() and self.spreading_rate_path!=None:
                synth = psk.make_synthetic(self.min_age,self.max_age,layer_depth,layer_thickness,layer_mag,azi,rd,ri,ad,ai,fix_sta,fix_end,twf,self.timescale_path,spreading_rate_path=self.spreading_rate_path,sz_name=self.dsk_row["sz_name"],length=syn_length,buff=self.syn_buff)
                srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
            else:
                synth = psk.make_synthetic(self.min_age,self.max_age,layer_depth,layer_thickness,layer_mag,azi,rd,ri,ad,ai,fix_sta,fix_end,twf,self.timescale_path,spreading_rate=spreading_rate,length=syn_length,buff=self.syn_buff)
                srf = lambda x,y: spreading_rate
        except: #Be wary of Errors here this is done so you can plot just the data and not have a million errors but it can mask true behavior be ready to need to debug this
            synth = [[0],[0],0]
            srf = lambda x,y: 0

        #Update Readouts on synthetic
        self.samp_dis_box.SetValue("%.2f"%float(synth[2]))

        #Update Data
        try:
            infile = os.path.join(self.dsk_row["data_dir"],self.dsk_row["comp_name"])
            if not os.path.isfile(infile): self.user_warning("Data file %s could not be found"%infile)
            self.dsk_row["strike"] = (azi-90)%360
            self.dsk_row["phase_shift"] = phase_shift
            self.deskew_df.set_value(self.dsk_idx,"strike",(azi-90)%360)
            self.deskew_df.set_value(self.dsk_idx,"phase_shift",phase_shift)

            #Center Synthetic
            if self.max_age!=self.min_age: step = (self.max_age-self.min_age)/(syn_length-1)
            else: step = .01
            srf_sz = lambda x: step*srf(self.dsk_row["sz_name"],x)
            dis_anom_min = sum(map(srf_sz,np.arange(0,self.dsk_row["age_min"]+step,step)))
            dis_anom_max = sum(map(srf_sz,np.arange(0,self.dsk_row["age_max"]+step,step)))
            anom_width = abs(dis_anom_max-dis_anom_min)/2
            center_dis = dis_anom_max - anom_width
            neg_anom = (-dis_anom_min + anom_width) - center_dis
            dis_synth = np.array(synth[1]) - center_dis
        except AttributeError: dis_synth,anom_width=synth[1],0

        ylim,xlim = self.ax.get_ylim(),self.ax.get_xlim()
        self.ax.clear()

#        psk.remove_axis_lines_and_ticks(self.ax)

        if self.m_show_anoms.IsChecked() or self.m_show_major_anoms.IsChecked(): self.plot_anomaly_names(major_anomolies_only=self.m_show_major_anoms.IsChecked(),anom_width=anom_width)

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
            else: self.user_warning("Cannot show other componenet for track type: %s"%str(self.dsk_row["track_type"]))

        try:
            psk.plot_skewness_data(self.dsk_row,self.dsk_row["phase_shift"],self.ax,color='k',zorder=3,picker=True)
            self.ax.axvspan(-anom_width,anom_width, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)
            if self.min_age<=-self.dsk_row["age_min"]: self.ax.axvspan(neg_anom-anom_width,neg_anom+anom_width, ymin=0, ymax=1.0, zorder=0, alpha=.5,color='yellow',clip_on=False,lw=0)
            self.ax.annotate("%s\n%s\n"%(self.dsk_row["sz_name"],self.track)+r"%.1f$^\circ$N,%.1f$^\circ$E"%(float(self.dsk_row['inter_lat']),utl.convert_to_0_360(self.dsk_row['inter_lon'])),xy=(0.02,1-0.15),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize)
        except AttributeError: pass

        self.ax.plot(dis_synth,synth[0],'r-',alpha=.4,zorder=1)
        self.ax.plot(dis_synth,np.zeros(len(dis_synth)),'k--')

#        psk.plot_scale_bars(self.ax,offset_of_bars = .05)

        scale = np.sqrt(sum(np.array(xlim)**2))
        if not scale<20 or scale>3000:
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
        self.track_box.SetItems([""]+self.deskew_df["comp_name"].tolist())
        self.on_select_track(-1)

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

    def on_open_ts_button(self,event):
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

    def on_open_dsk_button(self,event):
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

    def on_open_sr_file(self,event):
        dlg = wx.FileDialog(
            self, message="Choose SR Model File",
            defaultDir=self.WD,
            wildcard="Files (*.txt)|*.txt|All Files (*.*)|*.*",
            style=wx.FD_OPEN
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.spreading_rate_path=dlg.GetPath()
            self.sr_path.SetValue(self.spreading_rate_path)
            srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
            try: self.sr_box.SetValue("%.1f"%((srf(self.dsk_row["sz_name"],self.dsk_row["age_min"])+srf(self.dsk_row["sz_name"],self.dsk_row["age_max"]))/2))
            except AttributeError: pass
            if not self.m_use_sr_model.IsChecked(): self.m_use_sr_model.Check()
            try: self.srmw.update()
            except AttributeError: pass
        dlg.Destroy()

    def on_open_as_file(self,event):
        dlg = wx.FileDialog(
            self, message="Choose AS Model File",
            defaultDir=self.WD,
            wildcard="Files (*.txt)|*.txt|All Files (*.*)|*.*",
            style=wx.FD_OPEN
            )
        if dlg.ShowModal() == wx.ID_OK:
            self.anomalous_skewness_path=dlg.GetPath()
            self.as_path.SetValue(self.anomalous_skewness_path)
            if not self.m_use_as_model.IsChecked(): self.m_use_as_model.Check()
        dlg.Destroy()


    ##########################Drop Down Boxes################################

    def on_select_track(self,event):
        self.track = self.track_box.GetValue()
        if self.track=="" or self.track==None or self.track=="None": return
        self.dsk_row = self.deskew_df[self.deskew_df["comp_name"]==self.track].iloc[0]
        self.dsk_idx = self.deskew_df.index[self.deskew_df["comp_name"]==self.track]
        try:
            if "strike" in self.dsk_row and not np.isnan(float(self.dsk_row["strike"])):
                self.azi_box.SetValue("%.1f"%(float(self.dsk_row["strike"])-90))
            if "phase_shift" in self.dsk_row and not np.isnan(float(self.dsk_row["phase_shift"])):
                self.phase_shift_box.SetValue("%.1f"%float(self.dsk_row["phase_shift"]))
        except TypeError: self.user_warning("Invalid Strike or Phase Shift in deskew file for %s"%self.track)
        if self.m_use_sr_model.IsChecked() and self.spreading_rate_path!=None:
            srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
            self.sr_box.SetValue("%.1f"%((srf(self.dsk_row["sz_name"],self.dsk_row["age_min"])+srf(self.dsk_row["sz_name"],self.dsk_row["age_max"]))/2))
        if self.srmw_open: self.srmw.sz_box.SetValue(self.dsk_row["sz_name"]); self.srmw.on_select_sz(event)
        if self.tvw_open: self.tvw.on_parent_select_track()

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
        try: tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        except AttributeError: self.user_warning("No timescale file defined")
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

    def on_show_adjacent(self,event):
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
        event.Skip()

    def on_move_mouse_plot(self,event):
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)

        self.annotate_point(pos,xy=(1-0.12,1-0.11),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize)

        self.plot_tracer_point(pos[0],linestyle='--',color='red',alpha=.5)
        try:
            if self.tvw_open:
                self.tvw.plot_tracer_point(self.dsk_row,pos[0],color="red",marker="o",s=10)
                self.tvw.canvas.draw()
        except AttributeError: pass

        self.canvas.draw()
        event.Skip()

    def on_select_dleft_click(self,event):
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax.transData.inverted().transform(pos)
        self.set_new_intercept(pos[0])
        self.update(event)
        if self.tvw_open: self.tvw.update()

    ##########################Menu Functions################################

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
            if not self.m_use_sr_model.IsChecked() or self.spreading_rate_path==None:
                sr_value = self.sr_box.GetValue()
                srf = lambda x,y: sr_value
            else: srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
            if not self.m_use_as_model.IsChecked() or self.anomalous_skewness_path==None: asf = lambda x: 0
            else: asf = sk.generate_anomalous_skewness_model(self.anomalous_skewness_path)
            sk.create_max_file(self.deskew_df,srf,asf,outfile=outfile)
        dlg.Destroy()

    def on_save_maxtab_file(self,event):
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
            os.remove(".tmp.deskew")
        dlg.Destroy()

    def on_show_anoms(self,event):
        self.update(event)

    def on_show_major_anoms(self,event):
        self.update(event)

    def on_use_sr_model(self,event):
        self.update(event)

    def on_use_as_model(self,event):
        self.update(event)

    def on_srm_edit(self,event):
        if not self.srmw_open:
            try: starting_sz = self.dsk_row["sz_name"]
            except (AttributeError,KeyError) as e: starting_sz = ""
            self.srmw = SRMWindow(self.spreading_rate_path,parent=self,starting_sz=starting_sz,fontsize=self.fontsize,dpi=self.dpi)
            self.srmw.Center()
            self.srmw.Show()
            self.srmw_open=True

    def on_tv(self,event):
        if not self.tvw_open:
            self.tvw = TVWindow(parent=self,dpi=self.dpi)
            self.tvw.Center()
            self.tvw.Show()
            self.tvw_open=True

    def on_open_debug(self,event):
        pdb.set_trace()

    ##########################Additional Plotting and Backend Functions################

    def plot_anomaly_names(self,major_anomolies_only=False,anom_width=None):
        tdf = pd.read_csv(self.timescale_path,sep='\t',header=0,index_col=0)
        if self.m_use_sr_model.IsChecked() and self.spreading_rate_path!=None:
            srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
        else:
            sr = float(self.sr_box.GetValue())
            srf = lambda x,y: sr
        step = (self.max_age-self.min_age)/(float(self.samp_n_box.GetValue())-1)
        try: sz_name,srf_sz = self.dsk_row["sz_name"],lambda x: step*srf(sz_name,x)
        except (AttributeError,KeyError) as e: srf_sz = lambda x: step*srf("",x)
        if anom_width==None:
            dis_anom_min = sum(map(srf_sz,np.arange(0,self.dsk_row["age_min"]+step,step)))
            dis_anom_max = sum(map(srf_sz,np.arange(0,self.dsk_row["age_max"]+step,step)))
            anom_width = abs(dis_anom_max-dis_anom_min)/2
        central_center = 0
        major_anomolies_composed_of_subchrons = []
        for j,(i,row) in enumerate(tdf.iterrows()):
            try: i,base,top = str(i),float(row["base"]),float(row["top"])
            except ValueError: self.user_warning("Non-Numeric Chron Age for %s in %s"%(str(i),self.timescale_path)); return
            if major_anomolies_only and i.count('.')!=0:
                major_chron = i.split('.')[0]
                if major_chron in major_anomolies_composed_of_subchrons: continue
                else:
                    idxs = tdf.index[tdf.index.str.contains(major_chron)]
                    base = max(tdf.loc[idxs]["base"])
                    top = min(tdf.loc[idxs]["top"])
                    i = major_chron
                    major_anomolies_composed_of_subchrons.append(major_chron)
            elif base>self.max_age or top<self.min_age: continue
            anom_dis = sum(map(srf_sz,np.arange(0,base+step,step)))
            self.ax.axvline(central_center+anom_dis,linestyle='--',color='blue',alpha=.5)
            self.ax.axvline(central_center-anom_dis,linestyle='--',color='blue',alpha=.5)
            if j==0: self.ax.annotate(i,xy=(central_center,0),va="bottom",ha="center",fontsize=self.fontsize)
            else:
                anom_center = (anom_dis+(sum(map(srf_sz,np.arange(0,top+step,step)))))/2
                self.ax.annotate(i,xy=(central_center+anom_center,0),va="bottom",ha="center",fontsize=self.fontsize)
                self.ax.annotate(i,xy=(central_center-anom_center,0),va="bottom",ha="center",fontsize=self.fontsize)

    def plot_tracer_point(self,x,**kwargs):
        try: self.point_on_track.remove()
        except (AttributeError,ValueError) as e: pass
        self.point_on_track = self.ax.axvline(x,**kwargs)

    def annotate_point(self,pos,**kwargs):
        try: self.point_annotation.remove()
        except (AttributeError,ValueError) as e: pass
        self.point_annotation = self.ax.annotate("x = %.2f\ny = %.2f"%(float(pos[0]),float(pos[1])),**kwargs)

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

    def set_new_intercept(self,dis,dist_e=.5):
        new_lon,new_lat,cor_dis = sk.get_lon_lat_from_plot_pick(self.dsk_row,dis,dist_e=dist_e)
        if self.user_warning("Are you sure you want to correct site location by %.2f to coordinates (%.2f,%.2f)"%(cor_dis,new_lat,new_lon)):
            self.deskew_df.set_value(self.dsk_idx,"inter_lat",new_lat)
            self.deskew_df.set_value(self.dsk_idx,"inter_lon",new_lon)
            self.dsk_row = self.deskew_df.loc[self.dsk_idx].iloc[0]




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
