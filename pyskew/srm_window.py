import wx, os, sys
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.utilities as utl
import wx.lib.mixins.listctrl  as  listmix
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

class SRMWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,spreading_rate_path,parent=None,fontsize=8,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="SynthMagGUI V0.0.1",style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi=dpi
        self.spreading_rate_path = spreading_rate_path

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

        self.update()

    def init_UI(self):
        spacing = 10

        #---------------------------------Make ListCtrl for SR---------------------------------------------------#

        self.logger = EditableListCtrl(self.panel, ID=wx.ID_ANY, size=(300,300),style=wx.LC_REPORT)
        self.logger.InsertColumn(0, 'End Age',width=150)
        self.logger.InsertColumn(1, 'Half Rate',width=150)
#        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_click_listctrl, self.logger)
#        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK,self.on_right_click_listctrl,self.logger)
#        self.Bind(wx.EVT_LIST_ITEM_SELECTED, self.on_select_measurement, self.logger)

        #-----------------------------Make DropDown Box and Update-----------------------------------------------#

        sz_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Choose Spreading Zone"), wx.VERTICAL)

        self.sz_box = wx.ComboBox(self.panel, id=wx.ID_ANY,size=(150, 25),choices=[], style=wx.CB_DROPDOWN|wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_COMBOBOX, self.on_select_sz,self.sz_box)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_enter_sz, self.sz_box)

        update_button = wx.Button(self.panel, id=wx.ID_ANY, label='Save and Update',size=(150,25))
        self.Bind(wx.EVT_BUTTON, self.on_update_button, update_button)

        sz_sizer.AddMany([(self.sz_box, 1, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.TOP|wx.EXPAND,spacing),
                          (update_button, 1, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.TOP|wx.BOTTOM|wx.EXPAND,spacing)])

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

        #----------------------------------Build UI and Fit--------------------------------------------------------#

        side_bar_sizer = wx.BoxSizer(wx.VERTICAL)
        side_bar_sizer.AddMany([(sz_sizer,1,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.TOP|wx.EXPAND,spacing),
                                (self.canvas,4,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,spacing)])

        outer_sizer = wx.BoxSizer(wx.HORIZONTAL)
        outer_sizer.AddMany([(self.logger,1,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND),
                             (side_bar_sizer,3,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND)])

        self.panel.SetSizerAndFit(outer_sizer)

    def create_menu(self):
        pass #TODO

    #########################Update UI Funcions#############################

    def update(self):
        if self.spreading_rate_path==None: return

        

    def on_close_main(self,event):
        self.Destroy()

    ###################Button and Dropdown Functions#########################

    def on_select_sz(self,event):
        self.sz = self.sz_box.GetValue()

    def on_enter_sz(self,event):
        sz_tmp = self.sz_box.GetValue()
        if sz_tmp in sz_box.GetItems():
            self.sz = sz_tmp
        else:
            if self.sr_dict!=None and parent.user_warning("Spreading Zone provided not in spreading rate model would you like to add to model?"):
                self.sz = sz_tmp
                self.sr_dict[self.sz] = [1e10,40]
            else: self.sz_box.SetValue("");self.sz = None

    def on_update_button(self,event):
        pass #TODO

    ###########################Figure Funcions###############################

    def on_middle_click_plot(self,event):
        if event.LeftIsDown() or event.ButtonDClick():
            return
        elif self.plot_setting == "Zoom":
            self.plot_setting = "Pan"
            self.toolbar.pan('off')
        elif self.plot_setting == "Pan":
            self.plot_setting = "Zoom"
            self.toolbar.zoom()

    ###########################Utility Methods###############################

    def open_sr_model_for_sz(self,sz_name):
        sr_model = utl.open_sr_model_file(self.spreading_rate_path)
        self.sz_sr_model = sr_model[sz_name]
        return self.sz_sr_model

##################################################SubClasses#######################################################
class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
    ''' TextEditMixin allows any column to be edited. '''

    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """Constructor"""
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)
