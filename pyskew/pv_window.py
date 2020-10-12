import wx, os, sys
import numpy as np
import numpy.ma as ma
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.plot_geographic as pgeo
import pyskew.utilities as utl
from scipy.signal import hilbert
from geographiclib.geodesic import Geodesic
import wx.lib.buttons as buttons
import wx.lib.mixins.listctrl as listmix
import matplotlib as mpl
import matplotlib.path as mpath
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

class PVWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,dpi=100,geoid=Geodesic.WGS84,fontsize=8):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Phase Viewer %s"%parent.__version__,style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi = dpi
        self.geoid = geoid
        self.fontsize = fontsize
        self.mouse_left_down = False

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()
        self.configure()

        self.update()

    def init_UI(self):
        spacing = 10

        #-----------------------------------------Make SideBar-----------------------------------------------------#

        side_bar_sizer = wx.BoxSizer(wx.VERTICAL)

        #####-------------------------------------Make FilterBox---------------------------------------------------#

        filter_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Filter Options"), wx.VERTICAL)
        highlow_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Highcut/Lowcut"), wx.HORIZONTAL)
        order_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Order"), wx.HORIZONTAL)
        text_input_sizer = wx.BoxSizer(wx.HORIZONTAL)

        null_filter = lambda x,y,z,fs,order: x
        self.filters = {"None":null_filter,"Butterworth Bandpass":sk.butter_bandpass_filter,"Butterworth Lowpass":sk.butter_lowpass_filter}
        list_filters = ["None","Butterworth Bandpass","Butterworth Lowpass"]
        self.filter_type_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(200, 25), value=list_filters[0], choices=list_filters, style=wx.CB_DROPDOWN|wx.TE_READONLY)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_filter,self.filter_type_box)

        self.lowcut_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))
        self.highcut_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))
        self.order_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))

        order_sizer.Add(self.order_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)
        highlow_sizer.AddMany([(self.lowcut_box, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                               (self.highcut_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])
        text_input_sizer.AddMany([(highlow_sizer, 2, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                                  (order_sizer, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])
        filter_sizer.AddMany([(self.filter_type_box, 1, wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                              (text_input_sizer, 1, wx.ALIGN_CENTER|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL, spacing)])

        #####-------------------------------------Make BoundsBox---------------------------------------------------#

        bounds_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Bounds"), wx.HORIZONTAL)
        bounds_diff_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.low_bound_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))
        self.high_bound_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))
        self.aero_diff_box = wx.TextCtrl(self.panel, id=wx.ID_ANY|wx.TE_CENTRE, size=(50,25))

        bounds_sizer.AddMany([(self.low_bound_box, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                              (self.high_bound_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])
        bounds_diff_sizer.AddMany([(bounds_sizer, 2, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                                   (self.aero_diff_box, 1, wx.ALIGN_RIGHT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])

        #####-------------------------------------Make UpdateBtn---------------------------------------------------#

        self.update_button = wx.Button(self.panel, id=wx.ID_ANY, label='Update',size=(200,50))
        self.Bind(wx.EVT_BUTTON, self.on_update_button, self.update_button)

        #-------------------------------------Make Figures---------------------------------------------------------#

        self.power_fig = Figure((3, 3), dpi=self.dpi)
        self.power_canvas = FigCanvas(self.panel, -1, self.power_fig)
        self.power_toolbar = NavigationToolbar(self.power_canvas)
#        psk.remove_axis_lines_and_ticks(self.ax)
        self.power_toolbar.Hide()
        self.power_plot_setting = "Zoom"
        self.power_toolbar.zoom()
        self.power_canvas.Bind(wx.EVT_MIDDLE_DOWN,self.on_middle_click_power_plot)

        self.fig = Figure((4, 3), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        self.toolbar = NavigationToolbar(self.canvas)
#        psk.remove_axis_lines_and_ticks(self.ax)
        self.toolbar.Hide()
        self.plot_setting = "Zoom"
        self.toolbar.zoom()
        self.canvas.Bind(wx.EVT_MIDDLE_DOWN,self.on_middle_click_plot)
        self.canvas.Bind(wx.EVT_MOTION,self.on_move_mouse_plot)
        self.canvas.Bind(wx.EVT_LEFT_DOWN, self.on_left_click_down)
        self.canvas.Bind(wx.EVT_LEFT_UP, self.on_left_click_up)
#        self.canvas.Bind(wx.EVT_RIGHT_DCLICK, self.on_select_dright_click)

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        side_bar_sizer.AddMany([(filter_sizer, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                                (bounds_diff_sizer, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                                (self.update_button, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing),
                                (self.power_canvas, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND|wx.ALL, spacing)])

        outer_sizer = wx.BoxSizer(wx.HORIZONTAL)
        outer_sizer.AddMany([#(grd_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND|wx.LEFT|wx.RIGHT,spacing),
                             (side_bar_sizer,1,wx.ALIGN_CENTER|wx.ALIGN_TOP|wx.EXPAND),
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

    def configure(self):
        self.lowcut_box.SetValue("%.2f"%0.01)
        self.highcut_box.SetValue("%.2f"%0.1)
        self.order_box.SetValue("%d"%3)
        self.low_bound_box.SetValue("%.2f"%(-50.0))
        self.high_bound_box.SetValue("%.2f"%50.0)
        self.aero_diff_box.SetValue("%.1f"%0.0)


    #########################Update UI Funcions#############################

    def update(self):

        self.fig.clear() #clear
        self.power_fig.clear()

        self.draw_figures() #draw

        self.canvas.draw() #rerender
        self.power_canvas.draw()

    def on_close_main(self,event):
        self.parent.pv_open=False
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

    def on_select_filter(self,event):
        self.update()

    def on_update_button(self,event):
        self.update()

    ###########################Figure Funcions###############################

    def on_middle_click_plot(self,event):
        if event.LeftIsDown() or event.ButtonDClick(): event.Skip(); return
        elif self.plot_setting == "Zoom":
            self.plot_setting = "Pan"
            self.toolbar.pan('off')
        elif self.plot_setting == "Pan":
            self.plot_setting = "Pick"
            self.toolbar.pan('off')
            myCursor= wx.StockCursor(wx.CURSOR_IBEAM)
            self.canvas.SetCursor(myCursor)
        elif self.plot_setting == "Pick":
            self.plot_setting = "Zoom"
            self.toolbar.zoom()
        event.Skip()

    def on_middle_click_power_plot(self,event):
        if event.LeftIsDown() or event.ButtonDClick(): event.Skip(); return
        elif self.power_plot_setting == "Zoom":
            self.power_plot_setting = "Pan"
            self.power_toolbar.pan('off')
        elif self.power_plot_setting == "Pan":
            self.power_plot_setting = "Zoom"
            self.power_toolbar.zoom()
        event.Skip()

    def on_move_mouse_plot(self,event):
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax0.transData.inverted().transform(pos)

#        self.annotate_point(pos,xy=(1-0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize,va='top',ha='right')

        self.plot_tracer_point(pos[0],linestyle='--',color='red',alpha=.5)

        self.canvas.draw()

        event.Skip()

    def on_left_click_down(self,event):
        if self.plot_setting!="Pick": event.Skip(); return
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax0.transData.inverted().transform(pos)
        self.tmp_new_bounds = [pos[0],None]
        self.mouse_left_down = True
        event.Skip()

    def on_left_click_up(self,event):
        if self.plot_setting!="Pick": event.Skip(); return
        pos=event.GetPosition()
        width, height = self.canvas.get_width_height()
        pos = [pos[0],height-pos[1]]
        pos = self.ax0.transData.inverted().transform(pos)
        self.tmp_new_bounds[1] = pos[0]
        self.mouse_left_down = False

        if abs(self.tmp_new_bounds[1]-self.tmp_new_bounds[0])>0:
            self.low_bound_box.SetValue("%.1f"%float(min(self.tmp_new_bounds)))
            self.high_bound_box.SetValue("%.1f"%float(max(self.tmp_new_bounds)))

        self.update()

        event.Skip()

    ##########################Additional Plotting and Backend Functions################

    def on_parent_select_track(self):
        self.update()

    def plot_tracer_point(self,x,**kwargs):
        try:
            for tracer in self.tracers:
                tracer.remove()
        except (AttributeError,ValueError) as e: pass
        self.tracers = []
        self.tracers.append(self.ax0.axvline(x,**kwargs))
        self.tracers.append(self.ax1.axvline(x,**kwargs))
        self.tracers.append(self.ax2.axvline(x,**kwargs))
        if self.mouse_left_down:
            self.tracers.append(self.ax0.axvspan(min([self.tmp_new_bounds[0],x]),max([self.tmp_new_bounds[0],x]),color="yellow",alpha=.4))
            self.tracers.append(self.ax1.axvspan(min([self.tmp_new_bounds[0],x]),max([self.tmp_new_bounds[0],x]),color="yellow",alpha=.4))
            self.tracers.append(self.ax2.axvspan(min([self.tmp_new_bounds[0],x]),max([self.tmp_new_bounds[0],x]),color="yellow",alpha=.4))

    def draw_figures(self):

        ####################################################Get Values
        dsk_row = self.parent.dsk_row
        track = self.parent.track
        ddis = float(self.parent.samp_dis_box.GetValue())
        if ddis==0: self.parent.user_warning("Synthetic is required for comparision of phase, start by initalilzing a synthetic"); return
        synth_dis = self.parent.dis_synth
        synth_mag = self.parent.synth

        filter_type = self.filter_type_box.GetValue()
        lowcut = float(self.lowcut_box.GetValue())
        highcut = float(self.highcut_box.GetValue())
        order = int(self.order_box.GetValue())

        left_bound = float(self.low_bound_box.GetValue())
        right_bound = float(self.high_bound_box.GetValue())
        aero_diff = float(self.aero_diff_box.GetValue())

        left_idx = np.argmin(np.abs(synth_dis-left_bound))
        right_idx = np.argmin(np.abs(synth_dis-right_bound))
        left_idx,right_idx = min([left_idx,right_idx]),max([left_idx,right_idx])

        bin_range,bin_num = (-180,180),120

        ###################################################Filter Data

        data_path = os.path.join(dsk_row["data_dir"],dsk_row["comp_name"])
        data_df = utl.open_mag_file(data_path)
        projected_distances = utl.calc_projected_distance(dsk_row['inter_lon'],dsk_row['inter_lat'],data_df['lon'].tolist(),data_df['lat'].tolist(),(180+dsk_row['strike'])%360)
        shifted_mag = sk.phase_shift_data(data_df["mag"],dsk_row["phase_shift"])
        if np.any(np.diff(projected_distances["dist"])<0): itshifted_mag = np.interp(-synth_dis,-projected_distances["dist"],shifted_mag)
        else: itshifted_mag = np.interp(synth_dis,projected_distances["dist"],shifted_mag)
        fitshifted_mag = self.filters[filter_type](itshifted_mag,lowcut,highcut,fs=1/ddis,order=order)

        ###################################################Actual Plotting

        outer = gridspec.GridSpec(4, 1)

        ###################################################Axis 0: Magnetic profiles
        self.ax0 = self.fig.add_subplot(outer[0])

        if self.parent.show_other_comp: #Handle Other Aeromag Component

            if dsk_row["track_type"]=="aero":
                if "Ed.lp" in track:
                    other_track = track.replace("Ed.lp","Vd.lp")
                    other_phase = dsk_row["phase_shift"]-90
                elif "Hd.lp" in track:
                    other_track = track.replace("Hd.lp","Vd.lp")
                    other_phase = dsk_row["phase_shift"]-90
                elif "Vd.lp" in track:
                    other_track = track.replace("Vd.lp","Ed.lp")
                    if other_track not in self.parent.deskew_df["comp_name"].tolist(): other_track = track.replace("Vd.lp","Hd.lp")
                    other_phase = dsk_row["phase_shift"]+90
                else: self.parent.user_warning("Improperly named component files should have either Ed.lp, Hd.lp, or Vd.lp got: %s"%track); return
                oth_row = self.parent.deskew_df[self.parent.deskew_df["comp_name"]==other_track].iloc[0]
                psk.plot_skewness_data(oth_row,other_phase,self.ax0,xlims=[None,None],color='darkgreen',zorder=2,picker=True,alpha=.7,return_objects=True,flip=True)

                oth_data_path = os.path.join(oth_row["data_dir"],oth_row["comp_name"])
                oth_data_df = utl.open_mag_file(oth_data_path)
                oth_shifted_mag = sk.phase_shift_data(oth_data_df["mag"],other_phase)
                if np.any(np.diff(projected_distances["dist"])<0): oth_itshifted_mag = np.interp(-synth_dis,-projected_distances["dist"],oth_shifted_mag)
                else: oth_itshifted_mag = np.interp(synth_dis,projected_distances["dist"],oth_data_df)
                oth_fitshifted_mag = self.filters[filter_type](oth_itshifted_mag,lowcut,highcut,fs=1/ddis,order=order)
                self.ax0.plot(synth_dis,oth_fitshifted_mag,color="#299C29",zorder=3,alpha=.6)

            else: pass

        psk.plot_skewness_data(dsk_row,dsk_row["phase_shift"],self.ax0,xlims=[None,None],zorder=3,picker=True,return_objects=True,flip=True)
        self.ax0.plot(synth_dis,fitshifted_mag,color="#7F7D7D",zorder=3,alpha=.6)
        self.ax0.plot(self.parent.dis_synth,self.parent.synth,'r-',alpha=.4,zorder=1)
        self.ax0.set_ylabel("Magnetic Profiles")
#        self.ax0.get_xaxis().set_ticklabels([])

        ###################################################Axis 1/2: Phase Angles and Differences

        self.ax1 = self.fig.add_subplot(outer[1], sharex=self.ax0)
        self.ax2 = self.fig.add_subplot(outer[2], sharex=self.ax0)

        ###################################################Calculate: Phase Differences
        trimmed_dis = synth_dis[left_idx:right_idx]
        trimmed_synth = synth_mag[left_idx:right_idx]
        trimmed_fitshifted_mag = fitshifted_mag[left_idx:right_idx]

        al_data = np.angle(hilbert(fitshifted_mag),deg=False)[left_idx:right_idx]
        al_synth = np.angle(hilbert(np.real(synth_mag)),deg=False)[left_idx:right_idx]

        data_synth_diff = phase_diff_func(al_synth,al_data)

        if self.parent.show_other_comp and dsk_row["track_type"]=="aero":
            trimmed_oth_fitshifted_mag = oth_fitshifted_mag[left_idx:right_idx]
            al_oth = np.angle(hilbert(oth_fitshifted_mag),deg=False)[left_idx:right_idx]

            oth_synth_diff = phase_diff_func(al_synth,al_oth)
            oth_data_diff = phase_diff_func(al_oth,al_data)

            if abs(aero_diff) > 0:
                idx = ma.array(np.abs(oth_data_diff)<aero_diff)

                self.ax1.plot((trimmed_dis[~idx]),(np.rad2deg(al_oth[~idx])),color="darkgreen",linestyle=":")

                self.ax2.plot((trimmed_dis[~idx]),(oth_synth_diff[~idx]),color="tab:pink",alpha=.8,linestyle=":")
                self.ax2.plot((trimmed_dis[~idx]),(oth_data_diff[~idx]),color="tab:grey",alpha=.8,linestyle=":")

                self.ax1.plot((trimmed_dis[~idx]),(np.rad2deg(al_data[~idx])),color="k",linestyle=":")
                self.ax1.plot((trimmed_dis[~idx]),(np.rad2deg(al_synth[~idx])),color="r",linestyle=":")

                self.ax2.plot((trimmed_dis[~idx]),(data_synth_diff[~idx]),color="tab:red",alpha=.8,linestyle=":")

#                import pdb; pdb.set_trace()
#                not_trimmed_dis = (trimmed_dis[~idx])
#                not_trimmed_dis[np.diff(~idx,prepend=[0])] = ma.masked
#                not_al_data = (al_data[~idx])
#                not_al_data[np.diff(~idx)] = ma.masked
#                not_al_synth = (al_synth[~idx])
#                not_al_synth[np.diff(~idx)] = ma.masked
#                not_al_oth = (al_oth[~idx])
#                not_al_oth[np.diff(~idx)] = ma.masked
#                not_data_synth_diff = (data_synth_diff[~idx])
#                not_data_synth_diff[np.diff(~idx)] = ma.masked
#                not_oth_synth_diff = (oth_synth_diff[~idx])
#                not_oth_synth_diff[np.diff(~idx)] = ma.masked
#                not_oth_data_diff = (oth_data_diff[~idx])
#                not_oth_data_diff[np.diff(~idx)] = ma.masked
                trimmed_dis = (trimmed_dis[idx])
                al_data = (al_data[idx])
                al_synth = (al_synth[idx])
                al_oth = (al_oth[idx])
                data_synth_diff = (data_synth_diff[idx])
                oth_synth_diff = (oth_synth_diff[idx])
                oth_data_diff = (oth_data_diff[idx])

            self.ax1.plot(trimmed_dis,np.rad2deg(al_oth),color="darkgreen")

            self.ax2.plot(trimmed_dis,oth_synth_diff,color="tab:pink",alpha=.8)
            self.ax2.plot(trimmed_dis,oth_data_diff,color="tab:grey",alpha=.8)

        self.ax1.plot(trimmed_dis,np.rad2deg(al_data),color="k")
        self.ax1.plot(trimmed_dis,np.rad2deg(al_synth),color="r")

        self.ax2.plot(trimmed_dis,data_synth_diff,color="tab:red",alpha=.8)
#        self.ax2.get_xaxis.set_ticklabels
        self.ax0.set_xlim(*self.parent.ax.get_xlim())
        self.ax0.set_ylim(*self.parent.ax.get_ylim())

        self.ax1.set_ylabel("Phase Angles")
        self.ax2.set_ylabel("Phase Differences")


        ###################################################Axis 2.1: Power Spectrum
#        inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2])#, hspace=0.)
#        self.ax2 = self.fig.add_subplot(inner[0])


        ###################################################Axis 2.2: Phase Statistics
        self.ax3 = self.fig.add_subplot(outer[3])

        if self.parent.show_other_comp and dsk_row["track_type"]=="aero":
            self.ax3.hist(oth_synth_diff,range=bin_range,bins=bin_num,color="tab:pink",alpha=.5,zorder=2)
            self.ax3.hist(oth_data_diff,range=bin_range,bins=bin_num,color="tab:grey",alpha=.5,zorder=1)

        self.ax3.hist(data_synth_diff,range=bin_range,bins=bin_num,color="tab:red",alpha=.5,zorder=3)
        self.ax3.axvline(np.median(data_synth_diff),color="k",alpha=.5,zorder=5,linestyle=":")
        self.ax3.axvline(np.mean(data_synth_diff),color="k",alpha=.5,zorder=5)
        self.ax3.axvspan(np.mean(data_synth_diff)-np.std(data_synth_diff),np.mean(data_synth_diff)+np.std(data_synth_diff),color="tab:grey",alpha=.3,zorder=0)

        self.ax3.annotate(r"$\theta_{mean}$ = $%.1f^\circ \pm %.1f^\circ$"%(np.mean(data_synth_diff),np.std(data_synth_diff)) + "\n" + r"$\theta_{median}$ = %.1f$^\circ$"%np.median(data_synth_diff),xy=(0.02,1-0.02),xycoords="axes fraction",bbox=dict(boxstyle="round", fc="w",alpha=.5),fontsize=self.fontsize,va='top',ha='left')

        self.ax3.set_ylabel(r"$\Delta \theta$ Count")

        self.fig.suptitle("%s\n%s\n"%(dsk_row["sz_name"],track))

        ###################################################Power Figure

        synth_freqs = np.fft.fftfreq(len(synth_dis[left_idx:right_idx]),ddis)
        tdata_freqs = np.fft.fftfreq(len(shifted_mag[left_idx:right_idx]),ddis)
        tshifted_freq = np.fft.fft(shifted_mag[left_idx:right_idx])
        fitshifted_freq = np.fft.fft(fitshifted_mag[left_idx:right_idx])
        tsynth_freq = np.fft.fft(synth_mag[left_idx:right_idx])

        self.power_ax = self.power_fig.add_subplot(111)

        self.power_ax.plot(tdata_freqs, np.abs(tshifted_freq), color="k")
        self.power_ax.plot(synth_freqs, np.abs(fitshifted_freq), color="#7F7D7D")
        self.power_ax.plot(synth_freqs, np.abs(tsynth_freq), color="r")

        self.power_ax.set_xlim(0.0,max(synth_freqs))

def phase_diff_func(a1,a2):
#    a = np.rad2deg((rconvert_360((a1-a2)/2)))
    a = np.rad2deg((a1-a2))
    a = np.where(a<-180, a+360, a)
    a_cor = np.where(a>180, a-360, a)
    return a_cor
