# S2SNet - Sequences to Star Networks (2012)
# ============================================
# Calculate descriptors of sequence graphs
# Authors:
#  Cristian R Munteanu:    muntisa@gmail.com
#  Humberto González-Díaz: gonzalezdiazh@yahoo.es

# display a graph image in a Scrolled Windows

import wx

class MyScrolledWindow(wx.Frame):
    def __init__(self, parent, id, title, sFig):
        wx.Frame.__init__(self, parent, id, title, size=(710, 760))
        sw = wx.ScrolledWindow(self)
        bmp = wx.Image(sFig,wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        wx.StaticBitmap(sw, -1, bmp)
        wx.StaticText(sw, -1, '* You can find PNG pictures, DOT files and the plotting executables inside \"dot\" folder', (5, 680))
        sw.SetScrollbars(20, 20, 55, 40)#(20,20,55,40)
        
