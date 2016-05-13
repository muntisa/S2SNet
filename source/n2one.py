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

class Filters2(wx.Dialog):
    def __init__(self, parent, id, title):
        wx.Dialog.__init__(self, parent, id, title,size=(330,270))
        
        self.SetFont(wx.Font(8, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        #self.SetCursor(wx.StockCursor(wx.CURSOR_HAND))
        self.SetBackgroundColour('#929BF3')        
        
        ## FILES
        wx.StaticBox(self, -1, 'Files', (5, 15), size=(310, 160))

        wx.StaticText(self, -1, 'N-char Seqs', (25, 50))
        self.iSeqs = wx.TextCtrl(self, -1, 'cN-Seqs.txt',(85,45) , (120,-1) ,  style=wx.TE_LEFT)

        wx.StaticText(self, -1, 'Code', (25, 80))
        self.Code = wx.TextCtrl(self, -1, 'cN1-Code.txt',(85,75) , (120,-1) ,  style=wx.TE_LEFT)        

        wx.StaticText(self, -1, '1-char Seqs', (25, 110))
        self.fSeqs = wx.TextCtrl(self, -1, 'c1-Seqs.txt',(85,105) , (120,-1) ,  style=wx.TE_LEFT)

        wx.StaticText(self, -1, 'Groups', (25, 140))
        self.Grs = wx.TextCtrl(self, -1, 'c1-groups.txt',(85,135) , (120,-1) ,  style=wx.TE_LEFT)   
        
        ## BUTTONS
        br1=wx.Button(self, 20, '...', (210,45),(30,-1))
        br2=wx.Button(self, 21, '...', (210,75),(30,-1))
        br3=wx.Button(self, 22, '...', (210,105),(30,-1))
        br3=wx.Button(self, 23, '...', (210,135),(30,-1))

        edit1=wx.Button(self, 30, 'Edit', (250,45),(50,-1))
        edit2=wx.Button(self, 31, 'Edit', (250,75),(50,-1))
        edit3=wx.Button(self, 32, 'Edit', (250,105),(50,-1))
        edit3=wx.Button(self, 33, 'Edit', (250,135),(50,-1))

        self.bSubmit=GenBitmapTextButton(self, 25, wx.Bitmap(orig_dir+'images/gtk-font.png'), 'N2One', (45, 190))
        self.bHelp=GenBitmapTextButton(self, 26, wx.Bitmap(orig_dir+'images/gtk-help.png'), 'Help', (130, 190))
        self.bQuit=GenBitmapTextButton(self, 27, wx.Bitmap(orig_dir+'images/gtk-quit.png'), 'Cancel', (215, 190))

        self.bSubmit.SetBezelWidth(2)
        self.bHelp.SetBezelWidth(2)
        self.bQuit.SetBezelWidth(2)

        wx.EVT_BUTTON(self, 20, self.OnBriSeqs)
        wx.EVT_BUTTON(self, 21, self.OnBrCode)
        wx.EVT_BUTTON(self, 22, self.OnBrfSeqs)
        wx.EVT_BUTTON(self, 23, self.OnBrGrs)

        wx.EVT_BUTTON(self, 25, self.OnSubmit)
        wx.EVT_BUTTON(self, 26, self.OnHelp)
        wx.EVT_BUTTON(self, 27, self.OnQuit)
    
        wx.EVT_BUTTON(self, 30, self.OnEditiSeqs)
        wx.EVT_BUTTON(self, 31, self.OnEditCode)
        wx.EVT_BUTTON(self, 32, self.OnEditfSeqs)
        wx.EVT_BUTTON(self, 33, self.OnEditGrs)

        self.orig_dir = globals()["orig_dir"]

    def OnQuit(self, event):
        self.Close()
        
    def OnBriSeqs(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.iSeqs.SetValue(path)
        
    def OnBrCode(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.Code.SetValue(path)

    def OnBrfSeqs(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.fSeqs.SetValue(path)

    def OnBrGrs(self, event):
        dir = os.getcwd()
        open_dlg = wx.FileDialog(self, message='Choose a file', defaultDir=dir, defaultFile='', 
			style=wx.OPEN|wx.CHANGE_DIR)
        if open_dlg.ShowModal() == wx.ID_OK:
            path = open_dlg.GetPath()
            self.Grs.SetValue(path)

    def OnEditiSeqs(self, event):
        sequences=self.iSeqs.GetValue()
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

    def OnEditCode(self, event):
        sequences=self.Code.GetValue()
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

    def OnEditfSeqs(self, event):
        groups=self.fSeqs.GetValue()
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

    def OnEditGrs(self, event):
        weights=self.Grs.GetValue()
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

    def OnSubmit(self, event):
        # verify input file
        ListF=[]
        ListF.append(self.iSeqs.GetValue())
        ListF.append(self.Code.GetValue())
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
        print "N to 1-Character Sequence @ S2SNet 3.0 beta"
        print "-----------------------------------------------------"
        print '\nProcessing '+self.iSeqs.GetValue()+' ... ',
        N2One(self.iSeqs.GetValue(),self.Code.GetValue(),self.fSeqs.GetValue(),self.Grs.GetValue())
        print 'Done!\n'
        self.bSubmit.Enable(True)

def N2One(niSeqs,nCode,nfSeqs,nGrs):
    # read the code file
    fnCode = open(nCode,"r")
    lCodes = fnCode.readlines()
    fnCode.close()

    #create Group file
    fnGrs = open(nGrs,"w")

    chN=[]
    ch1=[]
    for lCode in lCodes : # for each Seq make all the calculations
        if lCode[-1]=="\n" : lCode = lCode[:-1] # removing strange characters
        if lCode[-1]=="\r" : lCode = lCode[:-1]
        Codes = lCode.split("\t") # code, chain, seq
        chN.append(Codes[0])
        ch1.append(Codes[1])
        fnGrs.write(str(Codes[1])+'\n')
    fnGrs.close()
    
    # read iSeqs
    fniSeqs = open(niSeqs,"r")
    liSeqss = fniSeqs.readlines()
    fniSeqs.close()
    
    # create fSeqs
    fnfSeqs = open(nfSeqs,"w")

    L1=''
    L2=''
    iSeq=''
    l=0
    for nline in liSeqss : # for each Seq make all the calculations
        if nline[-1]=="\n" : nline = nline[:-1] # removing strange characters
        if nline[-1]=="\r" : nline = nline[:-1]
        (L1,L2,iSeq) = nline.split("\t") # L1, L2, iSeqs
        fSeq=str(L1)+'\t'+str(L2)+'\t'
        scale=len(iSeq)/len(chN[0])
        dChs=len(chN[0])
        for i in range(scale):
            iS=iSeq[i*dChs:i*dChs+dChs]
            for f in range(len(chN)):
                fS=chN[f]
                if iS==fS:
                    fSeq+=ch1[f]
        if l != (len(liSeqss)-1):
            fSeq+='\n'
        fnfSeqs.write(fSeq)
        l+=1
    fnfSeqs.close()
    return
    

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
        vNo=zeros((lenData),dtype=float32) ## vNo=zeros(lenData,Float)
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
