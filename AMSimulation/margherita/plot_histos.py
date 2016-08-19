#!/usr/bin/env python

from ROOT import TFile, TH1F, TH2F, TProfile, TF1, TText, gPad, gROOT, gStyle
import numpy as np
import os

imgdir = "figures/"
histos = {}
figures = []

def make_figures(infile='histos.root'):
    if not os.path.exists(imgdir):
        os.makedirs(imgdir)

    tfile = TFile.Open(infile)
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

    gPad.DrawFrame(-1.,-1.,1.,1.)
    text_empty = TText(0., 0., "Empty")
    text_empty.SetTextAlign(22); text_empty.Draw()
    gPad.Print(imgdir+"empty.png")


    #print figures
    figures[:] = [
        'npc0',
        'npc1',
        'npc2',
        'npc3',
        'npc4',
        'npc5',
        'npc6',
        'npc7',
        'npc8',
        'npc9',
        'npc10',
        'npc11',
        'npc0_vs_par0_yfit2',
        'npc0_vs_par1_yfit2',
        'npc0_vs_par2_yfit2',
        'npc0_vs_par3_yfit2',
        'npc1_vs_par0_yfit2',
        'npc1_vs_par1_yfit2',
        'npc1_vs_par2_yfit2',
        'npc1_vs_par3_yfit2',
        'npc2_vs_par0_yfit2',
        'npc2_vs_par1_yfit2',
        'npc2_vs_par2_yfit2',
        'npc2_vs_par3_yfit2',
        'npc3_vs_par0_yfit2',
        'npc3_vs_par1_yfit2',
        'npc3_vs_par2_yfit2',
        'npc3_vs_par3_yfit2',
        'npc4_vs_par0_yfit2',
        'npc4_vs_par1_yfit2',
        'npc4_vs_par2_yfit2',
        'npc4_vs_par3_yfit2',
        'npc5_vs_par0_yfit2',
        'npc5_vs_par1_yfit2',
        'npc5_vs_par2_yfit2',
        'npc5_vs_par3_yfit2',
        'npc6_vs_par0_yfit2',
        'npc6_vs_par1_yfit2',
        'npc6_vs_par2_yfit2',
        'npc6_vs_par3_yfit2',
        'npc7_vs_par0_yfit2',
        'npc7_vs_par1_yfit2',
        'npc7_vs_par2_yfit2',
        'npc7_vs_par3_yfit2',
        'npc8_vs_par0_yfit2',
        'npc8_vs_par1_yfit2',
        'npc8_vs_par2_yfit2',
        'npc8_vs_par3_yfit2',
        'npc9_vs_par0_yfit2',
        'npc9_vs_par1_yfit2',
        'npc9_vs_par2_yfit2',
        'npc9_vs_par3_yfit2',
        'npc10_vs_par0_yfit2',
        'npc10_vs_par1_yfit2',
        'npc10_vs_par2_yfit2',
        'npc10_vs_par3_yfit2',
        'npc11_vs_par0_yfit2',
        'npc11_vs_par1_yfit2',
        'npc11_vs_par2_yfit2',
        'npc11_vs_par3_yfit2',

        'err0',
        'err1',
        'err2',
        'err3',
        'err0_vs_par0_yfit2',
        'err0_vs_par1_yfit2',
        'err0_vs_par2_yfit2',
        'err0_vs_par3_yfit2',
        'err1_vs_par0_yfit2',
        'err1_vs_par1_yfit2',
        'err1_vs_par2_yfit2',
        'err1_vs_par3_yfit2',
        'err2_vs_par0_yfit2',
        'err2_vs_par1_yfit2',
        'err2_vs_par2_yfit2',
        'err2_vs_par3_yfit2',
        'err3_vs_par0_yfit2',
        'err3_vs_par1_yfit2',
        'err3_vs_par2_yfit2',
        'err3_vs_par3_yfit2',

        'err0',
        'err1',
        'err2',
        'err3',
        'err0_vs_pt_yfit2',
        'err0_vs_eta_yfit2',
        'empty',
        'empty',
        'err1_vs_pt_yfit2',
        'err1_vs_eta_yfit2',
        'empty',
        'empty',
        'err2_vs_pt_yfit2',
        'err2_vs_eta_yfit2',
        'empty',
        'empty',
        'err3_vs_pt_yfit2',
        'err3_vs_eta_yfit2',
        'empty',
        'empty',

        'err4',
        'err5',
        'empty',
        'empty',
        'err4_vs_pt_yfit2',
        'err4_vs_eta_yfit2',
        'empty',
        'empty',
        'err5_vs_pt_yfit2',
        'err5_vs_eta_yfit2',
        'empty',
        'empty',
    ]
    return

def make_more_figures():
    def display_fit(h, xxmin, xxmax, fun="gaus"):
        status = 99
        if h.Integral() > 0 and (xxmax > xxmin and h.Integral(h.FindFixBin(xxmin),h.FindFixBin(xxmax))):
            status = int(h.Fit(fun, "Q", "", xxmin, xxmax))
        if status == 0:
            h.fit = h.GetFunction(fun)
        else:
            h.fit = TF1("fa1", fun + "(0)")
        h.fit.SetLineWidth(2); h.fit.SetLineColor(h.GetLineColor())

    def display_yfit(h, xxmin, xxmax, fun="gaus"):
        nbinsx, xmin, xmax = h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
        htitle, htitle2 = h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()
        h.yfit1 = h.ProjectionX(h.GetName()+"_yfit1"); h.yfit1.Reset()
        h.yfit2 = h.ProjectionX(h.GetName()+"_yfit2"); h.yfit2.Reset()
        h.yfit1.SetTitle(";"+htitle+";"+htitle2)
        h.yfit2.SetTitle(";"+htitle+";"+htitle2)

        for b in xrange(1,nbinsx+1):
            py = h.ProjectionY("_py", b, b)
            if py.Integral() > 50:
                display_fit(py, xxmin, xxmax)
                p1, p2 = py.fit.GetParameter(1), py.fit.GetParameter(2)
                e1, e2 = py.fit.GetParError(1), py.fit.GetParError(2)
                #is_fit_ok = (abs(p1 - py.GetMean())/py.GetRMS() < 1) and (1e-2 < p2/py.GetRMS() < 1e2)
                #if not is_fit_ok:
                #    p1, p2 = py.GetMean(), py.GetRMS()
                #    e1, e2 = py.GetMeanError(), py.GetRMSError()
                h.yfit1.SetBinContent(b, p1)
                h.yfit1.SetBinError(b, p2)
                h.yfit2.SetBinContent(b, p2)
                h.yfit2.SetBinError(b, e2)

    for hname, h in histos.iteritems():
        if h.ClassName() == "TH2F":
            ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
            yymin, yymax = h.GetMean(2)-4.0*h.GetRMS(2), h.GetMean(2)+4.0*h.GetRMS(2)
            display_yfit(h, yymin, yymax)

            for yfit in ["yfit1", "yfit2"]:
                h1 = getattr(h, yfit)
                if yfit == "yfit1":
                    h1.GetYaxis().SetTitle("%s" % h1.GetYaxis().GetTitle())
                    h1.SetMinimum(ymin); h1.SetMaximum(ymax)
                else:
                    h1.GetYaxis().SetTitle("#sigma(%s)" % h1.GetYaxis().GetTitle())
                    h1.SetMinimum(0.); h1.SetMaximum(0.1)
                h1.SetStats(0); h1.Draw()

                if hname.endswith("_vs_pt"):
                    gPad.SetLogx(True)
                else:
                    gPad.SetLogx(False)

                imgname = h1.GetName()
                gPad.Print(imgdir+imgname+".png")
    return

def serve_figures(htmlfile='index.html'):
    html = '''<html><head><link href="bootstrap.min.css" rel="stylesheet"></head><body><div class="container"><div class="row">%s</div></div></body></html>'''

    writeme = []
    #for fig in figures:
    #    s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig
    #    writeme.append(s)
    for fig in figures:
        s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig.replace("_yfit2", "_yfit1")
        writeme.append(s)
    for fig in figures:
        s = '''<div class="col-md-3"><img src="%s.png" class="img-responsive"></div>''' % fig.replace("_yfit2", "")
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

