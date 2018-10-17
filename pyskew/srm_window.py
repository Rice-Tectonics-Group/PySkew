import wx, os, sys
import numpy as np
import pyskew.skewness as sk
import pyskew.plot_skewness as psk
import pyskew.utilities as utl
import wx.lib.mixins.listctrl  as  listmix
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from functools import cmp_to_key

class SRMWindow(wx.Frame):

    #########################Init Funcions#############################

    def __init__(self,spreading_rate_path,parent=None,starting_sz="",fontsize=8,dpi=200):
        """Constructor"""
        #call init of super class
        default_style = wx.MINIMIZE_BOX | wx.MAXIMIZE_BOX | wx.RESIZE_BORDER | wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX | wx.CLIP_CHILDREN | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WS_EX_CONTEXTHELP | wx.FRAME_EX_CONTEXTHELP
        wx.Frame.__init__(self, parent, title="Spreading Rate Model V0.1.1",style=default_style, size=(400*2,300*2))
        self.Bind(wx.EVT_CLOSE, self.on_close_main)

        self.parent=parent
        self.dpi=dpi
        self.spreading_rate_path = spreading_rate_path

        self.panel = wx.Panel(self,-1,size=(400*2,300*2))

        #Populate UI and Menu
        self.init_UI()
        self.create_menu()

        self.sz=starting_sz
        try: self.open_srm()
        except: self.srm = None

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

        add_rate_sizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, "Add Rate to Model"), wx.HORIZONTAL)
        self.add_end_age_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25))
        self.add_half_rate_box = wx.TextCtrl(self.panel, id=wx.ID_ANY, size=(50,25))
        self.add_rate_btn = wx.Button(self.panel, id=wx.ID_ANY, label='Add Rate',size=(50,25))
        self.Bind(wx.EVT_BUTTON, self.on_add_rate_btn, self.add_rate_btn)
        add_rate_sizer.AddMany([(self.add_end_age_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.BOTTOM|wx.EXPAND,spacing),
                                (self.add_half_rate_box, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.BOTTOM|wx.EXPAND,spacing),
                                (self.add_rate_btn, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,spacing),])

        update_delete_sizer = wx.BoxSizer(wx.HORIZONTAL)
        delete_btn = wx.Button(self.panel, id=wx.ID_ANY, label='Delete Selected',size=(75,25))
        self.Bind(wx.EVT_BUTTON, self.on_delete_btn, delete_btn)

        update_button = wx.Button(self.panel, id=wx.ID_ANY, label='Save and Update',size=(75,25))
        self.Bind(wx.EVT_BUTTON, self.on_update_button, update_button)
        update_delete_sizer.AddMany([(update_button, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.EXPAND,spacing),
                                     (delete_btn, 1, wx.ALIGN_LEFT|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.EXPAND,spacing)])

        sz_sizer.AddMany([(self.sz_box, 2, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.TOP|wx.EXPAND,spacing),
                          (add_rate_sizer, 3, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.TOP|wx.EXPAND,spacing),
                          (update_delete_sizer, 2, wx.ALIGN_CENTER|wx.ALIGN_CENTER_VERTICAL|wx.LEFT|wx.RIGHT|wx.TOP|wx.BOTTOM|wx.EXPAND,spacing)])

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
                                (self.canvas,8,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.LEFT|wx.RIGHT|wx.BOTTOM|wx.EXPAND,spacing)])

        outer_sizer = wx.BoxSizer(wx.HORIZONTAL)
        outer_sizer.AddMany([(self.logger,1,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND),
                             (side_bar_sizer,3,wx.ALIGN_LEFT|wx.ALIGN_TOP|wx.EXPAND)])

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

    def update_self_and_parent(self):
        utl.write_sr_model_file(self.srm,self.spreading_rate_path)
        try:
            srf,_ = sk.generate_spreading_rate_model(self.spreading_rate_path)
            self.parent.sr_box.SetValue("%.1f"%((srf(self.parent.dsk_row["sz_name"],self.parent.dsk_row["age_min"])+srf(self.parent.dsk_row["sz_name"],self.parent.dsk_row["age_max"]))/2))
        except (AttributeError,KeyError) as e: pass

        self.update()
        self.parent.update(-1)

    def update(self): #Populates Logger and makes plot

        if self.parent.spreading_rate_path!=self.spreading_rate_path:
            self.spreading_rate_path = self.parent.spreading_rate_path
            self.open_srm()

        if self.spreading_rate_path==None: return

        self.logger.DeleteAllItems()
        points = [[],[]]
        if self.sz in self.srm.keys():
            prev_age = 0
            for i,(end_age,half_rate) in enumerate(self.srm[self.sz]):
                self.logger.InsertItem(i, "%.3f"%end_age)
                self.logger.SetItem(i, 1, "%.3f"%half_rate)

                if end_age>1e3:end_age=prev_age+10
                if i==0 and end_age>10: prev_age=end_age-10
                if prev_age==0 and end_age<=0: prev_age=end_age-10
                points[0].append(prev_age)
                points[0].append(end_age)
                points[1].append(half_rate)
                points[1].append(half_rate)
                prev_age = end_age

        self.ax.clear()

        self.ax.plot(points[0],points[1],'k-')

        self.canvas.draw()

    def on_close_main(self,event):
        self.parent.srmw_open=False
        self.Destroy()

    ###################Button and Dropdown Functions#########################

    def on_select_sz(self,event):
        self.sz = self.sz_box.GetValue()
        self.update()

    def on_enter_sz(self,event):
        sz_tmp = self.sz_box.GetValue()
        if sz_tmp in sz_box.GetItems():
            self.sz = sz_tmp
        else:
            if self.srm!=None and self.parent.user_warning("Spreading Zone provided not in spreading rate model would you like to add to model?"):
                self.sz = sz_tmp
                self.srm[self.sz] = [1e10,40]
            else: self.sz_box.SetValue("");self.sz = None
        self.update()

    def on_update_button(self,event): #Saves the edits and calls update on self and parent

        new_sz_srm = []
        for i in range(self.logger.GetItemCount()):
            try:
                new_sz_srm.append([float(self.logger.GetItemText(i,0)),float(self.logger.GetItemText(i,1))])
            except ValueError:
                self.parent.user_warning("Half Rate and Age to add must be numbers got %s,%s instead"%(self.logger.GetItemText(i,0),self.logger.GetItemText(i,1)))
                return
        self.srm[self.sz] = new_sz_srm

        self.update_self_and_parent()

    def on_add_rate_btn(self,event):

        try: new_end_age,new_half_rate = float(self.add_end_age_box.GetValue()),float(self.add_half_rate_box.GetValue())
        except ValueError: self.parent.user_warning("Half Rate and Age to add must be numbers got %s,%s instead"%(str(new_end_age),str(new_half_rate))); return

        if new_end_age in np.array(self.srm[self.sz])[:,0]:
            self.srm[self.sz] = [[new_end_age,new_half_rate] if self.srm[self.sz][i][0]==new_end_age else self.srm[self.sz][i] for i in range(len(self.srm[self.sz]))]
        else:
            self.srm[self.sz] += [[new_end_age,new_half_rate]]
        self.srm[self.sz].sort(key=cmp_to_key(lambda x,y: x[0]-y[0]))

        self.update_self_and_parent()

    def on_delete_btn(self,event):

        next_i,del_idxs = self.logger.GetNextSelected(-1),[]
        if next_i == -1: return
        while next_i != -1:
            del_idxs.append(next_i)
            next_i = self.logger.GetNextSelected(next_i)

        new_sz_srm = []
        for i in range(self.logger.GetItemCount()):
            if i in del_idxs: continue #Skip deleted data
            try:
                new_sz_srm.append([float(self.logger.GetItemText(i,0)),float(self.logger.GetItemText(i,1))])
            except ValueError:
                self.parent.user_warning("Half Rate and Age to add must be numbers got %s,%s instead"%(self.logger.GetItemText(i,0),self.logger.GetItemText(i,1)))
                return

        self.srm[self.sz] = new_sz_srm

        self.update_self_and_parent()

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

    ###########################Utility Funcions###############################

    def open_srm(self):
        self.srm = utl.open_sr_model_file(self.spreading_rate_path)
        self.sz_box.SetItems([""] + sorted(list(self.srm.keys())))
        self.sz_box.SetValue(self.sz)


##################################################SubClasses#######################################################
class EditableListCtrl(listmix.TextEditMixin,wx.ListCtrl):
    ''' TextEditMixin allows any column to be edited. '''

    #----------------------------------------------------------------------
    def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition,
                 size=wx.DefaultSize, style=0):
        """Constructor"""
        wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
        listmix.TextEditMixin.__init__(self)
