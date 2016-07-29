#!/usr/bin/env python

from ROOT import TFile, TH1F, TH2F, TProfile, gPad, gROOT, gStyle
import numpy as np
import os

imgdir = "figures/"
histos = {}
figures = []

def make_figures(infile='histos.root'):
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)

    tfile = TFile.Open("histos.root")
    tfile.cd()
    for h in tfile.GetListOfKeys():
        h = h.ReadObj()
        h.SetDirectory(0)
        histos[h.GetName()] = h
    tfile.Close()

    for hname, h in histos.iteritems():
        if h.ClassName() == "TH1F":
            h.Draw("hist")
            imgname = hname
            gPad.Print(imgdir+imgname+".png")
        elif h.ClassName() == "TH2F":
            h.Draw("COLZ")

            if hname.endswith("_vs_pt"):
                gPad.SetLogx(True)
            else:
                gPad.SetLogx(False)

            imgname = hname
            gPad.Print(imgdir+imgname+".png")

    #print figures
    figures[:] = [
        'err0',
        'err1',
        'err2',
        'err3',
        'err0_vs_par0_gaus2',
        'err1_vs_par0_gaus2',
        'err2_vs_par0_gaus2',
        'err3_vs_par0_gaus2',
        'err0_vs_par1_gaus2',
        'err1_vs_par1_gaus2',
        'err2_vs_par1_gaus2',
        'err3_vs_par1_gaus2',
        'err0_vs_par2_gaus2',
        'err1_vs_par2_gaus2',
        'err2_vs_par2_gaus2',
        'err3_vs_par2_gaus2',
        'err0_vs_par3_gaus2',
        'err1_vs_par3_gaus2',
        'err2_vs_par3_gaus2',
        'err3_vs_par3_gaus2',
        'npc4',
        'npc5',
        'npc11',
        'npc10',
        'npc4_vs_par0_gaus2',
        'npc5_vs_par0_gaus2',
        'npc11_vs_par0_gaus2',
        'npc10_vs_par0_gaus2',
        'npc4_vs_par1_gaus2',
        'npc5_vs_par1_gaus2',
        'npc11_vs_par1_gaus2',
        'npc10_vs_par1_gaus2',
        'npc4_vs_par2_gaus2',
        'npc5_vs_par2_gaus2',
        'npc11_vs_par2_gaus2',
        'npc10_vs_par2_gaus2',
        'npc4_vs_par3_gaus2',
        'npc5_vs_par3_gaus2',
        'npc11_vs_par3_gaus2',
        'npc10_vs_par3_gaus2',
        'npc0',
        'npc1',
        'npc2',
        'npc3',
        'npc0_vs_par0_gaus2',
        'npc1_vs_par0_gaus2',
        'npc2_vs_par0_gaus2',
        'npc3_vs_par0_gaus2',
        'npc0_vs_par1_gaus2',
        'npc1_vs_par1_gaus2',
        'npc2_vs_par1_gaus2',
        'npc3_vs_par1_gaus2',
        'npc0_vs_par2_gaus2',
        'npc1_vs_par2_gaus2',
        'npc2_vs_par2_gaus2',
        'npc3_vs_par2_gaus2',
        'npc0_vs_par3_gaus2',
        'npc1_vs_par3_gaus2',
        'npc2_vs_par3_gaus2',
        'npc3_vs_par3_gaus2',
        'npc6',
        'npc7',
        'npc8',
        'npc9',
        'npc6_vs_par0_gaus2',
        'npc7_vs_par0_gaus2',
        'npc8_vs_par0_gaus2',
        'npc9_vs_par0_gaus2',
        'npc6_vs_par1_gaus2',
        'npc7_vs_par1_gaus2',
        'npc8_vs_par1_gaus2',
        'npc9_vs_par1_gaus2',
        'npc6_vs_par2_gaus2',
        'npc7_vs_par2_gaus2',
        'npc8_vs_par2_gaus2',
        'npc9_vs_par2_gaus2',
        'npc6_vs_par3_gaus2',
        'npc7_vs_par3_gaus2',
        'npc8_vs_par3_gaus2',
        'npc9_vs_par3_gaus2',
        'err4',
        'err4',
        'err4_vs_pt_gaus2',
        'err4_vs_eta_gaus2',
        'err0',
        'err1',
        'err2',
        'err3',
        'err0_vs_pt_gaus2',
        'err1_vs_pt_gaus2',
        'err2_vs_pt_gaus2',
        'err3_vs_pt_gaus2',
        'err0_vs_eta_gaus2',
        'err1_vs_eta_gaus2',
        'err2_vs_eta_gaus2',
        'err3_vs_eta_gaus2',
    ]
    return

def make_more_figures():
    myptbins = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 800.0, 1000.0, 2000.0, 7000.0]

    for hname, h in histos.iteritems():
        if h.ClassName() == "TH2F":
            nbinsx, xmin, xmax = h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
            nbinsy, ymin, ymax = h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
            htitle2, htitle = h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()

            gaus1 = TH1F(hname+"_gaus1", ";"+htitle2+";"+htitle, nbinsx, xmin, xmax)
            gaus2 = TH1F(hname+"_gaus2", ";"+htitle2+";"+htitle, nbinsx, xmin, xmax)
            gaus1.SetMinimum(ymin/2.); gaus1.SetMaximum(ymax/2.)
            gaus2.SetMinimum(ymin/2.); gaus2.SetMaximum(ymax/2.)

            if hname.endswith("_vs_pt"):
                gaus1.SetBins(len(myptbins)-1, np.array(myptbins))
                gaus2.SetBins(len(myptbins)-1, np.array(myptbins))
                gaus1.logx = True; gaus1.logy = False
                gaus2.logx = True; gaus2.logy = False
            else:
                gaus1.logx = False; gaus1.logy = False
                gaus2.logx = False; gaus2.logy = False

            for b in xrange(1,nbinsx+1):
                py = h.ProjectionY("_px", b, b)
                if py.Integral() > 10:
                    py.Fit("gaus", "q")
                    fun = py.GetFunction("gaus")
                    p1, p2 = fun.GetParameter(1), fun.GetParameter(2)
                    e1, e2 = fun.GetParError(1), fun.GetParError(2)
                    gaus1.SetBinContent(b, p1)
                    gaus1.SetBinError(b, p2)
                    gaus2.SetBinContent(b, 0)
                    gaus2.SetBinError(b, p2)

            gaus1.SetStats(0); gaus1.Draw("e")
            gPad.SetLogx(gaus1.logx); gPad.SetLogy(gaus1.logy)
            imgname = gaus1.GetName()
            gPad.Print(imgdir+imgname+".png")
            gaus2.SetStats(0); gaus2.Draw("e")
            gPad.SetLogx(gaus2.logx); gPad.SetLogy(gaus2.logy)
            imgname = gaus2.GetName()
            gPad.Print(imgdir+imgname+".png")
    return

def serve_figures(htmlfile='index.html'):
    html = '''<html><head><link href="bootstrap.min.css" rel="stylesheet"></head><body><div class="container"><div class="row">%s</div></div></body></html>'''

    writeme = []
#    for fig in figures:
#        s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig
#        writeme.append(s)
    for fig in figures:
        s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig.replace("_gaus2", "_gaus1")
        writeme.append(s)
    for fig in figures:
        s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig.replace("_gaus2", "")
        writeme.append(s)
    writeme = '\n'.join(writeme)
    writeme = html % writeme

    if not os.path.isfile(imgdir+"bootstrap.min.css"):
        import subprocess
        subprocess.check_call(['wget', 'https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css'])
        subprocess.check_call(['mv', 'bootstrap.min.css', imgdir])

    with open(imgdir+htmlfile, 'w') as f:
        f.write(writeme)
    print "Check out the webpage at file://%s/%s" % (os.getcwd(), imgdir+htmlfile)
    return

def main():
    gROOT.LoadMacro("tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle();")
    gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
    gROOT.SetBatch(True)
    gStyle.SetOptStat(111110)
    gStyle.SetNdivisions(505, "XY")
    gStyle.SetPalette(57)  # kBird
    gStyle.SetNumberContours(100)

    make_figures()
    make_more_figures()
    serve_figures()

# ______________________________________________________________________________
if __name__ == "__main__":

    main()

