# S2SNet - Sequences to Star Networks (2012)
# ============================================
# Calculate descriptors of sequence graphs
# Authors:
#  Cristian R Munteanu:    muntisa@gmail.com
#  Humberto González-Díaz: gonzalezdiazh@yahoo.es

import wx
import wx.html
import os

orig_dir=os.getcwd()+"\\"

class MyHtmlPanel(wx.Panel):
    """
    class MyHtmlPanel inherits wx.Panel and adds a button and HtmlWindow
    """
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.SetFont(wx.Font(16, wx.DEFAULT, wx.NORMAL, wx.NORMAL))
        
        self.SetBackgroundColour("white")
        self.html1 = wx.html.HtmlWindow(self, id, pos=(0,30), size=(700,500))

        path=globals()["orig_dir"]+'help.htm' # read a HTM/HTML file
        self.html1.LoadPage(path)
