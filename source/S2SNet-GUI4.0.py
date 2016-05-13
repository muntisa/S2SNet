# -*- coding: cp1252 -*-
# S2SNet - Sequences to Star Networks (2012)
# ============================================
# Calculate descriptors of sequence graphs
# Authors:
#  Cristian R Munteanu:    muntisa@gmail.com
#  Humberto González-Díaz: gonzalezdiazh@yahoo.es

import wx
from wx.lib.buttons import GenBitmapTextButton

import os
import random
from numpy import *

orig_dir=os.getcwd()+"\\"

import S2SGf   # S2SG function
import htlmDoc # html help support
import figures # viewer for Graph pictures
import numbers # filter for numbers
import n2one   # filter n-character sequence

class MyMenu(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(370, 550),style=wx.CAPTION | wx.SYSTEM_MENU | wx.MINIMIZE_BOX | wx.CLOSE_BOX)
        self.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        self.SetBackgroundColour('#80BFF2')

        # menu
        menubar = wx.MenuBar()
        file = wx.Menu()
        open = wx.MenuItem(file, 101, '&New\tCtrl+N', 'Create new text file')
        open.SetBitmap(wx.Bitmap('images/gtk-new.png'))
        file.AppendItem(open)
        file.AppendSeparator()
        quit = wx.MenuItem(file, 105, '&Quit\tCtrl+Q', 'Quit S2SNet')
        quit.SetBitmap(wx.Bitmap('images/gtk-quit.png'))
        file.AppendItem(quit)
        
        calc = wx.Menu()    
        S2SG = wx.MenuItem(calc, 108, '&Sequence to Star Network\tCtrl+S', 'Transform Sequences to Star Network topological indices')
        S2SG.SetBitmap(wx.Bitmap('images/gtk-execute.png'))
        calc.AppendItem(S2SG)
        
        No2Seq = wx.MenuItem(calc, 107, '&Numbers to Sequence ...\tCtrl+N', 'Convert numbers in sequence')
        No2Seq.SetBitmap(wx.Bitmap('images/gtk-convert.png'))
        calc.AppendItem(No2Seq)

        N2One = wx.MenuItem(calc, 106, 'N to 1-&Character Sequence ...\tCtrl+C', 'Convert N-characters code in 1-character code')
        N2One.SetBitmap(wx.Bitmap('images/gtk-font.png'))
        calc.AppendItem(N2One)
        
        help = wx.Menu()
        helps = wx.MenuItem(help, 109, '&Help\tCtrl+H', 'How to use S2SNet')
        helps.SetBitmap(wx.Bitmap('images/gtk-help.png'))
        help.AppendItem(helps)
        
        about = wx.MenuItem(help, 110, '&About S2SNet\tCtrl+A', 'Informations about S2SNet and the authors')
        about.SetBitmap(wx.Bitmap('images/gtk-about.png'))
        help.AppendItem(about)
        
        menubar.Append(file, '&File')
        menubar.Append(calc, '&Calculations')
        menubar.Append(help, '&Help')
        
        self.SetMenuBar(menubar)
        self.statusbar = self.CreateStatusBar()        
        self.statusbar.SetStatusText("Transform your sequences in Star Network topological indices!")
        self.statusbar.SetBackgroundColour('RED')

        panel = wx.Panel(self, -1, (235, 20), (105, 130), style=wx.SUNKEN_BORDER)
        panel.SetBackgroundColour('#80BFF2')
        self.picture = wx.StaticBitmap(panel)
        self.picture.SetBitmap(wx.Bitmap(orig_dir+'images/logo2.bmp'))

        panel2 = wx.Panel(self, -1, (235, 175), (110, 20))
        panel2.SetBackgroundColour('#80BFF2')
        self.picture2 = wx.StaticBitmap(panel2)
        self.picture2.SetBitmap(wx.Bitmap(orig_dir+'images/wxPython.gif'))
        
        ## PARAMETERS group
        wx.StaticBox(self, -1, 'Parameters', (5, 15), size=(200, 175))
        
        self.emb=wx.CheckBox(self, -1 ,'Embedded', (25, 40))
        self.wts=wx.CheckBox(self, -1 ,'Weights', (125, 40))
        self.markov=wx.CheckBox(self, -1 ,'Markov normalization', (25, 70))
        
        wx.StaticText(self, -1, 'Matrix/Index power', (25, 100)) # 160
        self.power=wx.SpinCtrl(self, -1, '5', (125, 95), (50, -1), min=1, max=5)
        
        self.dot=wx.CheckBox(self, -1 ,'Network plots', (25, 130)) 
        self.det=wx.CheckBox(self, -1 ,'Output details', (25, 160)) #100
        
        self.markov.SetValue(True)
        
        ## FILES
        wx.StaticBox(self, -1, 'Files', (5, 210), size=(350, 205))
        wx.StaticText(self, -1, 'Sequences', (25, 235))
        self.sequences = wx.TextCtrl(self, -1, 'seqs.txt',(85,230) , (120,-1) ,  style=wx.TE_LEFT)
        wx.StaticText(self, -1, 'Groups', (25, 265))
        self.groups = wx.TextCtrl(self, -1, 'groups.txt',(85,260) , (120,-1) ,  style=wx.TE_LEFT)
        wx.StaticText(self, -1, 'Weights', (25, 295))
        self.weights = wx.TextCtrl(self, -1, 'weights.txt',(85,290) , (120,-1) ,  style=wx.TE_LEFT)
        wx.StaticText(self, -1, 'Results', (25, 325))
        self.results = wx.TextCtrl(self, -1, 'results.txt',(85,320) , (120,-1) ,  style=wx.TE_LEFT)
        wx.StaticText(self, -1, 'Details', (25, 355))
        self.details = wx.TextCtrl(self, -1, 'details.txt',(85,350) , (120,-1) ,  style=wx.TE_LEFT)

        SeqList=[]
        wx.StaticText(self, -1, 'Networks', (25, 385))
        self.Seqs = wx.ComboBox(self, -1, pos=(85,380), size=(100, -1), choices=SeqList, style=wx.CB_READONLY)
        self.Seqs.Enable(False)
        
        self.ExeList=['twopi','neato','dot','circo','fdp']
        self.EXE = wx.ComboBox(self, -1, pos=(190,380), size=(60, -1), choices=self.ExeList, style=wx.CB_READONLY)
        self.EXE.SetSelection(0)
        self.EXE.Enable(False)
        
        ## BUTTONS
        br1=wx.Button(self, 20, '...', (220,230), (30,-1))
        br2=wx.Button(self, 21, '...', (220,260), (30,-1))
        br3=wx.Button(self, 22, '...', (220,290), (30,-1))
        br4=wx.Button(self, 23, '...', (220,320), (30,-1))
        br5=wx.Button(self, 24, '...', (220,350), (30,-1))

        edit1=wx.Button(self, 30, 'Edit', (260,230))
        edit2=wx.Button(self, 31, 'Edit', (260,260))
        edit3=wx.Button(self, 32, 'Edit', (260,290))
        edit4=wx.Button(self, 33, 'Edit', (260,320))
        edit5=wx.Button(self, 34, 'Edit', (260,350))

        self.bSubmit=GenBitmapTextButton(self, 25, wx.Bitmap(orig_dir+'images/gtk-execute.png'), 'S2SNet', (15, 430))
        self.bHelp=GenBitmapTextButton(self, 26, wx.Bitmap(orig_dir+'images/gtk-help.png'), 'Help', (100, 430))
        self.bAbout=GenBitmapTextButton(self, 27, wx.Bitmap(orig_dir+'images/gtk-about.png'), 'About', (185, 430))
        self.bQuit=GenBitmapTextButton(self, 28, wx.Bitmap(orig_dir+'images/gtk-quit.png'), 'Quit', (270, 430))

        self.bSubmit.SetBezelWidth(2)
        self.bHelp.SetBezelWidth(2)
        self.bAbout.SetBezelWidth(2)
        self.bQuit.SetBezelWidth(2)

        self.bDisplay=wx.Button(self, 29, 'Display', (260,380))

        wx.EVT_BUTTON(self, 20, self.OnBrowseSeqs)
        wx.EVT_BUTTON(self, 21, self.OnBrowseGr)
        wx.EVT_BUTTON(self, 22, self.OnBrowseWts)
        wx.EVT_BUTTON(self, 23, self.OnBrowseRes)
        wx.EVT_BUTTON(self, 24, self.OnBrowseDet)

        wx.EVT_BUTTON(self, 25, self.OnSubmit)
        wx.EVT_BUTTON(self, 26, self.OnHelp)
        wx.EVT_BUTTON(self, 27, self.OnAbout)
        wx.EVT_BUTTON(self, 28, self.OnQuit)
        wx.EVT_BUTTON(self, 29, self.OnDisplay)

        wx.EVT_BUTTON(self, 30, self.OnEditSeqs)
        wx.EVT_BUTTON(self, 31, self.OnEditGr)
        wx.EVT_BUTTON(self, 32, self.OnEditWgs)
        wx.EVT_BUTTON(self, 33, self.OnEditRes)
        wx.EVT_BUTTON(self, 34, self.OnEditDet)

        self.Bind(wx.EVT_MENU, self.OnOpenEditor, id=101)
        self.Bind(wx.EVT_MENU, self.OnQuit, id=105)
        self.Bind(wx.EVT_MENU, self.OnSubmit, id=108)
        self.Bind(wx.EVT_MENU, self.OnNo2Seqs, id=107)
        self.Bind(wx.EVT_MENU, self.OnN2One, id=106)
        self.Bind(wx.EVT_MENU, self.OnHelp, id=109)
        self.Bind(wx.EVT_MENU, self.OnAbout, id=110)

        # header
        print "\n**********************************************************"
        print "S2SNet - Sequence to Star Network (ver. 4.0, 2012)\n\nby Cristian R. Munteanu and Humberto Gonzalez-Diaz"
        print "\nE-mail: muntisa@gmail.com"
        print "**********************************************************\n"
        self.SeqList=[]
        self.orig_dir = globals()["orig_dir"] #takes the original folder for the un modified files in the browse controls
        
        self.bDisplay.Enable(False)
        
    def OnBrowseSeqs(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.sequences.SetValue(path)

    def OnBrowseGr(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.groups.SetValue(path)

    def OnBrowseWts(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.weights.SetValue(path)

    def OnBrowseRes(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.results.SetValue(path)

    def OnBrowseDet(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.details.SetValue(path)
            
    def OnSubmit(self, event):
        emb=self.emb.GetValue()
        wts=self.wts.GetValue()
        dot=self.dot.GetValue() 
        markov=self.markov.GetValue()
        det=self.det.GetValue()
        power=self.power.GetValue()
        sequences=self.sequences.GetValue()
        groups=self.groups.GetValue()
        weights=self.weights.GetValue()
        results=self.results.GetValue()
        details=self.details.GetValue()

        if sequences.find("\\")==-1:
            sequences=self.orig_dir+sequences 
        if groups.find("\\")==-1:
            groups=self.orig_dir+groups
        if weights.find("\\")==-1:
            weights=self.orig_dir+weights
        if results.find("\\")==-1:
            results=self.orig_dir+results
        if details.find("\\")==-1:
            details=self.orig_dir+details

        ListF=[]
        ListF.append(sequences)
        ListF.append(groups)
        if wts == True:
            ListF.append(weights)
        er=0
        for sFile in ListF:
            try:
                fsFile = open(sFile,"r")
                fsFile.close()
            except IOError, error:
                dlg = wx.MessageDialog(self, 'Error opening file "'+sFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                er=1
            except UnicodeDecodeError, error:
                dlg = wx.MessageDialog(self, 'The file '+sFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                er=1
        if er == 1: #if there is any open file error, the function is finishing
            return

        self.Seqs.Enable(False)
        self.EXE.Enable(False)
        self.bDisplay.Enable(False)
        # self.bSubmit.Enable(False)
        self.statusbar.SetStatusText('Transform Sequences in Star Network indices ... ')
        # S2SGmain returns the list of Seqs for Combo control
        nErr=0
        try:
            (self.SeqList,statGrAll,g)=S2SGf.S2SGmain(emb, wts, dot, markov, det, power, sequences, groups, weights,results,details,self.orig_dir)
            lenStatGrAll=len(statGrAll)
        except:
            nErr=1
            print "Please check the input file and restart the calculation."
            self.statusbar.SetStatusText('S2SNet error. Please check the input files.')
            return

        if self.dot.GetValue()==True and nErr!=1 and lenStatGrAll!=0: # only if DOT is selected!
            print '\nCreate MIN, MAX and AVG Networks:\n'
            # calculates min, max and avg of all the graphs!
            minGr=statGrAll[0].copy()
            maxGr=statGrAll[0].copy()
            k=len(minGr)
            nSeq=len(statGrAll)
            avgGr=zeros((k),dtype=float32) ## avgGr=zeros([k],Float)
            sumGr=zeros((k),dtype=float32)
            for iGr in range(nSeq):
                curGr=statGrAll[iGr]
                for iiGr in range(k):
                    if curGr[iiGr] < minGr[iiGr]:
                        minGr[iiGr]=curGr[iiGr]
                    else:
                        if curGr[iiGr] > maxGr[iiGr]:
                            maxGr[iiGr]=curGr[iiGr]
                    sumGr[iiGr]+=curGr[iiGr]
            avgGr=sumGr/nSeq
                    
            # create min, max and avg Seqs in a new Seqs file
            # MIN Graph
            sAllGraphs=self.orig_dir+'seqs_MinMaxAvg.txt'
            sResGraphs=self.orig_dir+'results_MinMaxAvg.txt'
            fallFile = open(sAllGraphs,"w")
            # removed temporally due to the ZERO division error
##            statSeq=''
##            for ik in range(len(minGr)):
##                temp=g[ik]
##                for istat in range(int(minGr[ik])):
##                    statSeq+=str(temp[0]) # create fake Seq with the first character of each group
##            fallFile.write('minNet\tAll\t'+statSeq)
            fallFile.close()
            # MAX Graph
            fallFile = open(sAllGraphs,"a")
            statSeq=''
            for ik in range(len(maxGr)):
                temp=g[ik]
                for istat in range(int(maxGr[ik])):
                    statSeq+=str(temp[0]) # create fake Seq with the first character of each group
            #fallFile.write('\nmaxNet\tAll\t'+statSeq)
            fallFile.write('maxNet\tAll\t'+statSeq)
            fallFile.close()
            # AVG Graph
            fallFile = open(sAllGraphs,"a")
            statSeq=''
            for ik in range(len(minGr)):
                temp=g[ik]
                for istat in range(int(avgGr[ik])):
                    statSeq+=str(temp[0]) # create fake Seq with the first character of each group
            fallFile.write('\navgNet\tAll\t'+statSeq)
            fallFile.close()
            # recalculate the min, max, avg indices and DOT files
            (self.SeqList2,statGrAll2,g2)=S2SGf.S2SGmain(False, False, True, markov, False, power, sAllGraphs, groups, weights,sResGraphs,details,self.orig_dir)

            self.Seqs.Clear() # clear combo for only calculations values
            self.Seqs.AppendItems(self.SeqList)
            self.Seqs.AppendItems(self.SeqList2)
            self.Seqs.SetSelection(0) # initialize with the first Seq from the list
            self.Seqs.Enable(True)
            self.EXE.Enable(True)
            self.bDisplay.Enable(True)
        self.statusbar.SetStatusText('S2SNet calculation finished')
        #self.bSubmit.Enable(True)
        
    def OnDisplay(self, event): # display PNGs
        self.bDisplay.Enable(False)
        wDir=orig_dir+"dot\\"
        item=self.Seqs.GetSelection()
        if item > len(self.SeqList)-1:
            item=item-len(self.SeqList)
            sFig=self.SeqList2[item]+'.png'
        else:
            sFig=self.SeqList[item]+'.png'
        sEXE=self.ExeList[self.EXE.GetSelection()]
        
        sType='png'
        sDot=sFig[:-4]+'.txt'
        
        self.statusbar.SetStatusText('Creating '+sFig+' ... ')

        S2SGf.Dot2Figure(sEXE, sType, sDot, sFig, wDir) #create pictures

        x,y = random.randint(1,50), random.randint(1,50)
        try:
            fDisplay = figures.MyScrolledWindow(None, -1, 'S2SNet Viewer: '+sFig+' - made with '+sEXE,wDir+sFig)
            fDisplay.SetIcon(wx.Icon(self.orig_dir+'favicon.ico', wx.BITMAP_TYPE_ICO))
            fDisplay.SetPosition((x,y))
            self.statusbar.SetStatusText('Image '+sFig+' was created in folder <dot>')
            fDisplay.Show(True)
        except:
            dlg = wx.MessageDialog(self, 'Error creating file "'+sFig+'". Please check if you have all the DLLs and EXEs inside DOT folder and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            self.statusbar.SetStatusText(sFig+' error')
            
        self.bDisplay.Enable(True)

    def OnNo2Seqs(self, event):
        dNumbers = numbers.Filters(None, -1, 'Numbers to Sequence @ S2SNet 4.0')
        dNumbers.SetIcon(wx.Icon(self.orig_dir+'favicon.ico', wx.BITMAP_TYPE_ICO))

        dS=wx.DisplaySize()
        dD=dNumbers.GetSize()
        pos=(dS[0]-dD[0]-20,20)
        dNumbers.SetPosition(pos)
        
        val = dNumbers.ShowModal()
        dNumbers.Destroy()        

    def OnN2One(self, event):
        dN2One = n2one.Filters2(None, -1, 'N to 1-Character Sequence @ S2SNet 4.0')
        dN2One.SetIcon(wx.Icon(self.orig_dir+'favicon.ico', wx.BITMAP_TYPE_ICO))

        dS=wx.DisplaySize()
        dD=dN2One.GetSize()
        pos=(dS[0]-dD[0]-20,dS[1]-dD[1]-50)
        dN2One.SetPosition(pos)
        
        val = dN2One.ShowModal()
        dN2One.Destroy()
        
    def OnOpenEditor(self, event):
        cmd = 'notepad.exe' # open  notepad
        os.system(cmd)

    def OnEditSeqs(self, event):
        sequences=self.sequences.GetValue()
        if sequences.find("\\")==-1:
            sequences=self.orig_dir+sequences 
        iFile=sequences        
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnEditGr(self, event):
        groups=self.groups.GetValue()
        if groups.find("\\")==-1:
            groups=self.orig_dir+groups 
        iFile=groups
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnEditWgs(self, event):
        weights=self.weights.GetValue()
        if weights.find("\\")==-1:
            weights=self.orig_dir+weights 
        iFile=weights
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnEditRes(self, event):
        results=self.results.GetValue()
        if results.find("\\")==-1:
            results=self.orig_dir+results
        iFile=results
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()

    def OnEditDet(self, event):
        details=self.details.GetValue()
        if details.find("\\")==-1:
            details=self.orig_dir+details
        iFile=details
        try:
            cmd = 'notepad.exe '+iFile # open  notepad
            os.system(cmd)
        except IOError, error:
            dlg = wx.MessageDialog(self, 'Error opening file "'+iFile+'". Please check if the file exists and try again.','S2SNet: File Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
        except UnicodeDecodeError, error:
            dlg = wx.MessageDialog(self, 'The file '+iFile+' contains at least one non-unicode characther. Please check your file and try again.' + str(error), 'S2SNet: Unicode Error', wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            
    def OnHelp(self, event):
        frame3 = wx.Frame(None, -1, "S2SNet Help", size=(720, 560),style=wx.CAPTION | wx.SYSTEM_MENU | wx.MINIMIZE_BOX | wx.CLOSE_BOX)
        frame3.SetIcon(wx.Icon(self.orig_dir+'favicon.ico', wx.BITMAP_TYPE_ICO))
        frame3.SetPosition((5,5))
        htlmDoc.MyHtmlPanel(frame3,-1)
        frame3.Show(True)
        
    def OnAbout(self, event):
        description = """
S2SNet is an application that is tranforming any sequence of characthers in Star Network topological indices:
- Shannon Entropy of the Markov Matrices;
- Connectivity matrix Traces;
- Harary number;
- Wiener index;
- Gutman topological index;
- Schultz topological index (non-trivial part);
- Moreau-Broto indices;
- Balaban distance connectivity index;
- Randic connectivity index and Kier-Hall indices.

It is a free Python application, with wxPython GUI and Graphviz graphics back-end, initially used for bioinformatics protein QSAR analysis.
"""

        licence = """
Free software for academinc use
S2SNet ver. 4.0
Copyright © S2SNet 2012
"""
        info = wx.AboutDialogInfo()
        info.SetIcon(wx.Icon(self.orig_dir+'images/logo.png', wx.BITMAP_TYPE_PNG))
        info.SetName('S2SNet - Sequence to Star Network')
        info.SetVersion('4.0')
        info.SetDescription(description)
        info.SetCopyright('© 2012 S2SNet')
        ## info.SetWebSite('http://miaja.tic.udc.es/software/S2SNet/')
        info.SetLicence(licence)
        info.AddDeveloper('Cristian R. Munteanu (muntisa@gmail.com)')
        info.AddDeveloper('\nHumberto Gonzalez-Diaz (gonzalezdiazh@yahoo.es)')
        wx.AboutBox(info)
        
    def OnQuit(self, event):
        self.Close()

class MyApp(wx.App):
    def OnInit(self):
        frame = MyMenu(None, -1, 'S2SNet-Sequence to Star Network ver.4.0')
        frame.SetIcon(wx.Icon('favicon.ico', wx.BITMAP_TYPE_ICO))
        frame.Centre()
        frame.Show(True)
        return True

app = MyApp(0)
app.MainLoop()
