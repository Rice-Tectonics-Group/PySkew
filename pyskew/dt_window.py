import wx, os, sys
import numpy as np
import pyskew.utilities as utl
import pyskew.skewness as sk

class DetrendWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,parent=None,fontsize=8,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Detrend Tool V0.1.1",style=default_style, size=(200,200))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)
        self.deg,self.pols,self.poly,self.projected_distances = 0,None,None,None

        self.parent=parent
        self.dpi=dpi

        self.panel = wx.Panel(self,-1,size=(200,200))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

    def init_UI(self):
        spacing = 10

        #-----------------------------Make TextCtrls for Fit Info-----------------------------------------------#

        deg_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Degree of Fit"), wx.HORIZONTAL)
        self.deg_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25))
        self.deg_box.SetValue(str(self.deg))
#        self.add_half_rate_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25))
        deg_sizer.AddMany([(self.deg_box, 1, wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.LEFT|wx.BOTTOM,spacing)])

        #---------------------------Make Buttons for Plot and Update--------------------------------------------#

        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)

        plot_btn = wx.Button(self.panel, id=wx.ID_ANY, label='Plot Fit',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_plot_btn, plot_btn)

        sub_btn = wx.Button(self.panel, id=wx.ID_ANY, label='Subtract Fit',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_sub_btn, sub_btn)

        btn_sizer.AddMany([(plot_btn, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.EXPAND,spacing),
                           (sub_btn, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.EXPAND,spacing)])

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        outer_sizer = wx.BoxSizer(wx.VERTICAL)
        outer_sizer.AddMany([(deg_sizer,1,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.ALL,spacing),
                             (btn_sizer,2,wx.ALIGN_LEFT|wx.ALIGN_BOTTOM|wx.EXPAND|wx.ALL,spacing)])

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

    def on_close_main(self,event):
        self.parent.dtw_open=False
        self.parent.update()
        self.Destroy()

    ###################Button Functions#########################

    def on_plot_btn(self,event):
        self.fit_poly()
        try: self.line.remove()
        except AttributeError as e: pass
        [self.line] = self.parent.ax.plot(self.projected_distances,self.poly,color="tab:orange")
        self.parent.canvas.draw()

    def on_sub_btn(self,event):
        if self.parent.user_warning("WARNING: This is a non-reversable function please be sure this is what you want to do. Also if you are not looking at the 0 phase profile this may have unexpected results, please view the 0 phase profile before hitting OK."):
            try: self.deg = int(self.deg_box.GetValue())
            except (ValueError,TypeError) as e: self.parent.user_warning("Degree of polynomial must be a natural number")
            try:
                mag_path = os.path.join(self.parent.dsk_row["data_dir"],self.parent.dsk_row["comp_name"])
                mag_df = utl.open_mag_file(mag_path)
                self.projected_distances = utl.calc_projected_distance(self.parent.dsk_row["inter_lon"],self.parent.dsk_row["inter_lat"],mag_df["lon"],mag_df["lat"],self.parent.dsk_row["strike"])["dist"]
                shifted_mag = mag_df['mag'].tolist()
                self.pols = np.polyfit(self.projected_distances,shifted_mag,self.deg)
                self.poly = np.polyval(self.pols,self.projected_distances)
                mag_df["mag"] = mag_df["mag"] - self.poly
                utl.write_mag_file_df(mag_df,mag_path)
            except AttributeError: return
        self.on_plot_btn(event)

    ###########################Utility Funcions###############################

    def fit_poly(self):
        try: self.deg = int(self.deg_box.GetValue())
        except (ValueError,TypeError) as e: self.parent.user_warning("Degree of polynomial must be a natural number")
        try:
            mag_path = os.path.join(self.parent.dsk_row["data_dir"],self.parent.dsk_row["comp_name"])
            mag_df = utl.open_mag_file(mag_path)
            self.projected_distances = utl.calc_projected_distance(self.parent.dsk_row["inter_lon"],self.parent.dsk_row["inter_lat"],mag_df["lon"],mag_df["lat"],self.parent.dsk_row["strike"])["dist"]
            shifted_mag = sk.phase_shift_data(mag_df['mag'].tolist(),self.parent.dsk_row["phase_shift"])
            self.pols = np.polyfit(self.projected_distances,shifted_mag,self.deg)
            self.poly = np.polyval(self.pols,self.projected_distances)
        except AttributeError: return


