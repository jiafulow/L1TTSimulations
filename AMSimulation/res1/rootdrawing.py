import sys
tmpargv = sys.argv[:]  # [:] for a copy, not reference
sys.argv = []

import ROOT
from ROOT import TH1, TH1F, TH2F, TProfile, TProfile2D, TFile, TFileCollection, TChain, TCanvas, TPad, TLegend, TLatex, TLine, TPolyLine, TBox, TGraph, TGraphErrors, TGraphAsymmErrors, TEfficiency, gROOT, gSystem, gStyle, gPad
from rootcolors import *
sys.argv = tmpargv

import re
import os
from math import sqrt, pi, ceil, floor
from random import randint
from itertools import izip
from array import array
import argparse
import tempfile
import json

# ______________________________________________________________________________
# Classes

class MyDrawer:
    def __init__(self):
        # ROOT
        gROOT.LoadMacro("tdrstyle.C")
        gROOT.ProcessLine("setTDRStyle()")

        gStyle.SetEndErrorSize(2)
        gStyle.SetPadRightMargin(0.05)
        gStyle.SetTitleOffset(1.1, "Y")
        gStyle.SetLabelSize(0.04, "XY")
        gStyle.SetLabelSize(0.04, "Z")
        gStyle.SetNdivisions(505, "XY")

        #gStyle.SetPalette(55)  # rainbow color map
        gStyle.SetPalette(kBird)  # new default color map
        gStyle.SetNumberContours(100)

        gStyle.SetOptStat(111110)
        gStyle.SetStatX(0.94)
        gStyle.SetStatY(0.93)
        gStyle.SetStatH(0.30)
        gStyle.SetStatW(0.28)

        TH1.SetDefaultSumw2()

class MyParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        outdir = (self.parser.prog.replace("drawer_", "figures_"))[:-3]
        self.parser.add_argument("infile", help="input file (either .root file or .txt file listing .root files)")
        self.parser.add_argument("--outdir", default=outdir, help="output directory (default: %(default)s)")
        self.parser.add_argument("--logx", action="store_true", help="draw plots with log x scale (default: %(default)s)")
        self.parser.add_argument("--logy", action="store_true", help="draw plots with log y scale (default: %(default)s)")
        self.parser.add_argument("-n", "--nentries", type=int, default=-1, help="number of entries (default: %(default)s)")
        self.parser.add_argument("-b", "--batch", action="store_true", help="batch mode without graphics (default: %(default)s)")
        self.parser.add_argument("-v", "--verbose", action="count", default=0, help="verbosity (default: %(default)s)")

    def add_argument(self, *args, **kwargs):
        self.parser.add_argument(*args, **kwargs)

    def set_defaults(self, **kwargs):
        self.parser.set_defaults(**kwargs)

    def parse_args(self, args=None, namespace=None):
        self.options = self.parser.parse_args(args=args, namespace=namespace)

        # Create outdir if necessary
        if not self.options.outdir.endswith("/"):
            self.options.outdir += "/"
        if not os.path.exists(self.options.outdir):
            os.makedirs(self.options.outdir)

        # Make input file list
        if not self.options.infile.endswith(".root") and not self.options.infile.endswith(".txt"):
            raise ValueError("infile must be .root file or .txt file")

        if self.options.infile.endswith(".root"):
            # Create a temporary file with one line
            with tempfile.NamedTemporaryFile() as infile:
                infile.write(self.options.infile)
                infile.flush()
                self.options.tfilecoll = TFileCollection("fc", "", infile.name)
        else:
            self.options.tfilecoll = TFileCollection("fc", "", self.options.infile)

        # Batch mode
        if self.options.batch:
            gROOT.SetBatch(True)
        return self.options


# ______________________________________________________________________________
# Globals

tlatex = TLatex()
tlatex.SetNDC()
tlatex.SetTextFont(42)
tlatex.SetTextSize(0.03)

tlegend = TLegend(0.70,0.74,0.96,0.94)
tlegend.SetFillStyle(0)
tlegend.SetLineColor(0)
tlegend.SetShadowColor(0)
tlegend.SetBorderSize(0)

tline = TLine()
tline.SetLineColor(kGray+2)
tline.SetLineStyle(2)

col = kBlack
fcol = kRed

donotdelete = []  # persist in memory


# ______________________________________________________________________________
# Math

def deltaPhi(phi1, phi2):
    result = phi1 - phi2
    while result > pi:  result -= 2*pi
    while result <= -pi:  result += 2*pi
    return result

def deltaEta(eta1, eta2):
    return eta1 - eta2

def deltaR(eta1, phi1, eta2, phi2):
    deta = deltaEta(eta1, eta2)
    dphi = deltaPhi(phi1, phi2)
    return sqrt(deta*deta + dphi*dphi)


# ______________________________________________________________________________
# Auxiliary

def count_if(seq, condition):
    return sum(1 for item in seq if condition(item))

def modify_binning(nbins, xmin, xmax):
    binsize = (xmax - xmin) / nbins
    return (nbins + 1, xmin - (binsize/2.), xmax + (binsize/2.))

def moveLegend(tlegend, x1, y1, x2, y2):
    # SetX1NDC, SetX2NDC etc don't update if the tlegend is not drawn first
    # SetX1, SetX2 etc don't update if the tlegend has already been drawn
    # Thanks ROOT!
    tlegend.SetX1(x1); tlegend.SetY1(y1); tlegend.SetX2(x2); tlegend.SetY2(y2)
    tlegend.SetX1NDC(x1); tlegend.SetY1NDC(y1); tlegend.SetX2NDC(x2); tlegend.SetY2NDC(y2)

def movePalette(h, x1=0.88, y1=0.13, x2=0.92, y2=0.95):
    paletteObj = h.FindObject("palette")
    paletteObj.SetX1NDC(x1); paletteObj.SetY1NDC(y1); paletteObj.SetX2NDC(x2); paletteObj.SetY2NDC(y2)

def moveStats(h, x1=0.64, y1=0.70, x2=0.88, y2=0.93):
    statsObj = h.FindObject("stats")
    statsObj.SetX1NDC(x1); statsObj.SetY1NDC(y1); statsObj.SetX2NDC(x2); statsObj.SetY2NDC(y2)

def CMS_label():
    cmsTextFont = 62; extraTextFont = 52
    lumiTextSize = 0.6; lumiTextOffset = 0.2; extraOverCmsTextSize = 0.76
    cmsTextSize = 0.75; cmsTextOffset = 0.1; extraTextSize = extraOverCmsTextSize * cmsTextSize
    relPosX = 0.045; relPosY = 0.035; relExtraDY = 1.2
    lumiText = "";  cmsText = "CMS";  extraText = "Preliminary Simulation"

    old = (tlatex.GetTextFont(), tlatex.GetTextSize())
    l, r, t, b = gPad.GetLeftMargin(), gPad.GetRightMargin(), gPad.GetTopMargin(), gPad.GetBottomMargin()
    posX_ = l + relPosX*(1-l-r); posY_ = 1-t+lumiTextOffset*t
    relPosX = 0.105; posX_ = 0.7 - relPosX*(1-l-r);  # right aligned
    tlatex.SetTextFont(42); tlatex.SetTextSize(lumiTextSize*t)
    tlatex.DrawLatex(1-r,1-t+lumiTextOffset*t, lumiText)
    tlatex.SetTextFont(cmsTextFont); tlatex.SetTextSize(cmsTextSize*t)
    tlatex.DrawLatex(posX_,1-t+lumiTextOffset*t,cmsText)  # out of frame
    tlatex.SetTextFont(extraTextFont); tlatex.SetTextSize(extraTextSize*t)
    tlatex.DrawLatex(posX_ + relPosX*(1-l-r), posY_, extraText)
    tlatex.SetTextFont(old[0]); tlatex.SetTextSize(old[1])
    return

def CMS_0T_label():
    cmsTextFont = 62; extraTextFont = 52; extraText2Font = 42
    lumiTextSize = 0.6; lumiTextOffset = 0.2; extraOverCmsTextSize = 0.76
    cmsTextSize = 0.75; cmsTextOffset = 0.1; extraTextSize = extraOverCmsTextSize * cmsTextSize;
    relPosX = 0.045; relPosY = 0.035; relExtraDY = 1.2
    lumiText = "";  cmsText = "CMS";  extraText = "Preliminary Simulation"; extraText2 = ", B = 0T"

    old = (tlatex.GetTextFont(), tlatex.GetTextSize())
    l, r, t, b = gPad.GetLeftMargin(), gPad.GetRightMargin(), gPad.GetTopMargin(), gPad.GetBottomMargin()
    posX_ = l + relPosX*(1-l-r); posY_ = 1-t+lumiTextOffset*t
    relPosX = 0.105; posX_ = 0.7 - relPosX*(1-l-r) - 0.1;  # right aligned
    tlatex.SetTextFont(42); tlatex.SetTextSize(lumiTextSize*t)
    tlatex.DrawLatex(1-r,1-t+lumiTextOffset*t, lumiText)
    tlatex.SetTextFont(cmsTextFont); tlatex.SetTextSize(cmsTextSize*t)
    tlatex.DrawLatex(posX_,1-t+lumiTextOffset*t,cmsText)  # out of frame
    tlatex.SetTextFont(extraTextFont); tlatex.SetTextSize(extraTextSize*t)
    tlatex.DrawLatex(posX_ + relPosX*(1-l-r), posY_, extraText)
    tlatex.SetTextFont(extraText2Font); tlatex.SetTextSize(extraTextSize*t)
    tlatex.DrawLatex(posX_ + relPosX*(1-l-r) + 0.25, posY_, extraText2)
    tlatex.SetTextFont(old[0]); tlatex.SetTextSize(old[1])
    return

def getMaximum(histos):
    maxima = []
    for h in histos:
        maxima.append(h.h.GetMaximum())
    return max(maxima)

def save(imgdir, imgname, redraw_axis=True, dot_pdf=True, dot_root=False, dot_c=False):
    if redraw_axis:
        gPad.RedrawAxis()
    #gPad.Modified(); gPad.Update()
    gPad.Print(imgdir+imgname+".png")
    if dot_pdf:
        gPad.Print(imgdir+imgname+".pdf")
    if dot_root:
        def save_for_root(p):
            if p.ClassName() == "TPad" or p.ClassName() == "TCanvas":
                for pp in p.GetListOfPrimitives():
                    save_for_root(pp)
            p.Write()
            if hasattr(p, 'additional'):
                for pp in p.additional:
                    pp.Write()

        tfile = TFile.Open(imgdir+imgname+".root", "RECREATE")
        save_for_root(gPad)
        tfile.Close()
    if dot_c:
        gPad.Print(imgdir+imgname+".C")

# ______________________________________________________________________________
# Useful functions

convert_key_to_int = lambda pairs: dict([(int(k),v) for (k,v) in pairs])

def get_reverse_map(direct_map):
    reverse_map = {}
    for i in xrange(6*8):
        for m in direct_map[i]:
            reverse_map.setdefault(m, []).append(i)
    return reverse_map

