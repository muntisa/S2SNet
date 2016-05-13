# S2SNet - Sequences to Star Networks (2012)
# ============================================
# Calculate descriptors of sequence graphs
# Authors:
#  Cristian R Munteanu:    muntisa@gmail.com
#  Humberto González-Díaz: gonzalezdiazh@yahoo.es

from numpy import *
import os

import wx
from wx.lib.buttons import GenBitmapTextButton

import htlmDoc # html help support

orig_dir=os.getcwd()+"\\"

class Filters(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title,size=(340,400))
        
        self.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        #self.SetCursor(wx.StockCursor(wx.CURSOR_HAND))
        self.SetBackgroundColour('#AAC3EF')
        
        ## PARAMETERS group
        wx.StaticBox(self, -1, 'Parameters', (50, 15), size=(240, 120))

        wx.StaticText(self, -1, 'Minimum', (80, 45))
        self.MinV = wx.TextCtrl(self, -1, '0',(140,40) , (50,-1) ,  style=wx.TE_LEFT)

        wx.StaticText(self, -1, 'Maximum', (80, 75))
        self.MaxV = wx.TextCtrl(self, -1, '100',(140,70) , (50,-1) ,style=wx.TE_LEFT)
        
        wx.StaticText(self, -1, 'Groups No.', (80, 105))
        self.GrsNo=wx.SpinCtrl(self, -1, '10', (140, 100), (50, -1), min=2, max=80)
        wx.StaticText(self, -1, '(max. 80)', (200, 105))
        
        ## FILES
        wx.StaticBox(self, -1, 'Files', (5, 150), size=(320, 160))
        
        wx.StaticText(self, -1, 'Numbers', (25, 180))
        self.numFile = wx.TextCtrl(self, -1, 'numbers.txt',(85,175) , (120,-1) ,  style=wx.TE_LEFT)

        wx.StaticText(self, -1, 'Sequences', (25, 210))
        self.sequences = wx.TextCtrl(self, -1, 'nseqs.txt',(85,205) , (120,-1) ,  style=wx.TE_LEFT)        

        wx.StaticText(self, -1, 'Groups', (25, 240))
        self.groups = wx.TextCtrl(self, -1, 'ngroups.txt',(85,235) , (120,-1) ,  style=wx.TE_LEFT)

        wx.StaticText(self, -1, 'Intervals', (25, 270))
        self.codif = wx.TextCtrl(self, -1, 'nintervals.txt',(85,265) , (120,-1) ,  style=wx.TE_LEFT)
        
        ## BUTTONS
        br1=wx.Button(self, 20, '...', (210,175),(30,-1))
        br2=wx.Button(self, 21, '...', (210,205),(30,-1))
        br3=wx.Button(self, 22, '...', (210,235),(30,-1))
        br3=wx.Button(self, 23, '...', (210,265),(30,-1))

        edit1=wx.Button(self, 30, 'Edit', (250,175),(50,-1))
        edit2=wx.Button(self, 31, 'Edit', (250,205),(50,-1))
        edit3=wx.Button(self, 32, 'Edit', (250,235),(50,-1))
        edit3=wx.Button(self, 33, 'Edit', (250,265),(50,-1))

        self.bSubmit=GenBitmapTextButton(self, 25, wx.Bitmap(orig_dir+'images/gtk-convert.png'), 'No2Seqs', (45, 320))
        self.bHelp=GenBitmapTextButton(self, 26, wx.Bitmap(orig_dir+'images/gtk-help.png'), 'Help', (130, 320))
        self.bQuit=GenBitmapTextButton(self, 27, wx.Bitmap(orig_dir+'images/gtk-quit.png'), 'Cancel', (215, 320))

        self.bSubmit.SetBezelWidth(2)
        self.bHelp.SetBezelWidth(2)
        self.bQuit.SetBezelWidth(2)
        
        getStat=wx.Button(self, 34, 'G\nE\nT', (200,40),(15,50),wx.CENTER)

        wx.EVT_BUTTON(self, 20, self.OnBrowseNum)
        wx.EVT_BUTTON(self, 21, self.OnBrowseSeqs)
        wx.EVT_BUTTON(self, 22, self.OnBrowseGr)
        wx.EVT_BUTTON(self, 23, self.OnBrowseCod)

        wx.EVT_BUTTON(self, 25, self.OnSubmit)
        wx.EVT_BUTTON(self, 26, self.OnHelp)
        wx.EVT_BUTTON(self, 27, self.OnQuit)
    
        wx.EVT_BUTTON(self, 30, self.OnEditNum)
        wx.EVT_BUTTON(self, 31, self.OnEditSeqs)
        wx.EVT_BUTTON(self, 32, self.OnEditGr)
        wx.EVT_BUTTON(self, 33, self.OnEditCod)
        wx.EVT_BUTTON(self, 34, self.OnGetStats)

        self.orig_dir = globals()["orig_dir"]

    def OnQuit(self, event):
        self.Close()
        
    def OnBrowseNum(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.numFile.SetValue(path)
        
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

    def OnBrowseCod(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.codif.SetValue(path)

    def OnEditNum(self, event):
        sequences=self.numFile.GetValue()
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

    def OnEditCod(self, event):
        weights=self.codif.GetValue()
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
            
    def OnHelp(self, event):
        frameH = wx.Frame(None, -1, "S2SNet Help", size=(720, 560),style=wx.CAPTION | wx.SYSTEM_MENU | wx.MINIMIZE_BOX | wx.CLOSE_BOX)
        frameH.SetIcon(wx.Icon(self.orig_dir+'favicon.ico', wx.BITMAP_TYPE_ICO))
        frameH.SetPosition((5,5))
        htlmDoc.MyHtmlPanel(frameH,-1)
        frameH.Show(True)

    def OnGetStats(self, event):
        # verify input file
        ListF=[]
        ListF.append(self.numFile.GetValue())
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
        #end verification
        (nMin,nMax,nAvg)=No2SeqStats(self.numFile.GetValue())
        self.MinV.SetValue(str(nMin))
        self.MinV.Refresh()
        self.MaxV.SetValue(str(nMax))
        self.MaxV.Refresh()

    def OnSubmit(self, event):
        # verify input file
        ListF=[]
        ListF.append(self.numFile.GetValue())
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
        #end verification
        self.bSubmit.Enable(False)
        # header
        #
        print "\n-----------------------------------------------------"
        print "Numbers to Sequence @ S2SNet 3.0 beta"
        print "-----------------------------------------------------"
        print '\nProcessing '+self.numFile.GetValue()+' ... '
        No2Seq(self.MinV.GetValue(),self.MaxV.GetValue(),self.GrsNo.GetValue(),self.numFile.GetValue(),self.sequences.GetValue(),self.groups.GetValue(),self.codif.GetValue())
        print '\nDone!\n'
        self.bSubmit.Enable(True)

def No2Seq(nMin,nMax,nGrs,nFile,seqFile,grsFile,codeFile): # min, max and no. of groups
    chars='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890()=*+-_:;.,><[]{}' # 80 ch

    # read numbers file
    fnFile = open(nFile,"r")
    nlines = fnFile.readlines()
    fnFile.close()
    
    SeqsNo=len(nlines)
    
    # calculate the groups
    niGrs=int(nGrs)
    U=abs(float(nMax)-float(nMin))/niGrs #the unit of one group

    # create groups and codification files
    fCodFile = open(codeFile,"w")
    fGrsFile = open(grsFile,"w")
    
    vFrom=zeros((niGrs)) ## vFrom=zeros(niGrs,Float) # has groups no items
    vTo=zeros((niGrs)) # has groups no items
    vCh={}
    for ig in range(niGrs):
        vFrom[ig]=(ig)*U+float(nMin)
        vTo[ig]=(ig+1)*U+float(nMin)
        vCh[ig]=chars[ig]
        fCodFile.write(vCh[ig]+'\t'+str(vFrom[ig])+'\t'+str(vTo[ig])+'\n')
        fGrsFile.write(vCh[ig]+'\n')
    fCodFile.close()
    fGrsFile.close()
    
    l=0
    Mins=zeros((SeqsNo)) ## Mins=zeros([SeqsNo],Float)
    Maxs=zeros((SeqsNo))
    Avgs=zeros((SeqsNo))

    # create Sequence file
    fSeqFile = open(seqFile,"w")
    
    for nline in nlines : # for each Seq make all the calculations
        if nline[-1]=="\n" : nline = nline[:-1] # removing strange characters
        if nline[-1]=="\r" : nline = nline[:-1]
        spLine = nline.split("\t") # code, chain, seq
        L1=spLine[0]
        L2=spLine[1]
        print '\t'+str(L1)+'\t'+str(L2)
        lenData=len(spLine[2:])
        Seq=''
        vNo=zeros((lenData))
        for i in range(lenData):
            vNo[i]=float(spLine[i+2])
            for ig in range(niGrs):
                if ig<(niGrs-1):
                    if vFrom[ig] <= vNo[i] and vNo[i] < vTo[ig]:
                        Seq=Seq+vCh[ig]
                if ig == (niGrs-1):
                    if ((vFrom[ig] <= vNo[i]) and (vNo[i] < vTo[ig])) or (vNo[i] == vTo[ig]):
                        Seq=Seq+vCh[ig]
        if l < SeqsNo-1:
            fSeqFile.write(str(L1)+'\t'+str(L2)+'\t'+Seq+'\n')
        else:
            fSeqFile.write(str(L1)+'\t'+str(L2)+'\t'+Seq)
        
        l+=1
    
    fSeqFile.close()
    return

def No2SeqStats(nFile):    
    # read numbers file
    fnFile = open(nFile,"r")
    nlines = fnFile.readlines()
    fnFile.close()
    
    l=0
    SeqsNo=len(nlines)
    Mins=zeros((SeqsNo))
    Maxs=zeros((SeqsNo))
    Avgs=zeros((SeqsNo))
    
    for nline in nlines : # for each Seq make all the calculations
        if nline[-1]=="\n" : nline = nline[:-1] # removing strange characters
        if nline[-1]=="\r" : nline = nline[:-1]
        spLine = nline.split("\t") # code, chain, seq
        L1=spLine[0]
        L2=spLine[1]
        lenData=len(spLine[2:])
        vNo=zeros((lenData))
        for i in range(lenData):
            vNo[i]=float(spLine[i+2])
        Mins[l]=min(vNo)
        Maxs[l]=max(vNo)
        Avgs[l]=float(average(vNo))
        l+=1
    
    return ( min(Mins), max(Maxs),average(Avgs) )
