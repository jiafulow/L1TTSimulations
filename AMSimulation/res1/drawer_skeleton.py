#!/usr/bin/env python

from rootdrawing import *


# ______________________________________________________________________________
# Drawer
def drawer_book(options):
    histos = {}

    nbinsx, xmin, xmax = 50, -1, 1
    hname = "h"
    histos[hname] = TH1F(hname, "; x", nbinsx, xmin, xmax)

    # Style
    for hname, h in histos.iteritems():
        h.style = 0; h.logx = options.logx; h.logy = options.logy
    donotdelete.append(histos)
    return histos

# ______________________________________________________________________________
def drawer_project(tree, histos, options):
    #tree.SetBranchStatus("*", 0)

    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break
        if (ievt % 1000 == 0):  print "Processing event: %i" % ievt

        hname = "h"
        histos[hname].Fill(1)

    #tree.SetBranchStatus("*", 1)
    return

# ______________________________________________________________________________
def drawer_draw(histos, options):
    for hname, h in histos.iteritems():

        if h.style == 0:  # filled
            h.SetLineWidth(2); h.SetMarkerSize(0)
            h.SetLineColor(col); h.SetFillColor(fcol)
        elif h.style == 1:  # marker
            h.SetLineWidth(2); h.SetMarkerStyle(20); h.SetFillStyle(0)
            h.SetLineColor(col); h.SetMarkerColor(col);
        if h.logy:
            h.SetMaximum(h.GetMaximum() * 14); h.SetMinimum(0.5)
        else:
            h.SetMaximum(h.GetMaximum() * 1.4); h.SetMinimum(0.)

        h.Draw("hist")
        gPad.SetLogy(h.logy)
        CMS_label()
        save(options.outdir, hname)
    return

# ______________________________________________________________________________
def drawer_sitrep(histos, options):
    print "--- SITREP --------------------------------------------------------"
    print


# ______________________________________________________________________________
# Main function
def main(options):

    # Init
    drawer = MyDrawer()
    tchain = TChain("ntupler/tree", "")
    tchain.AddFileInfoList(options.tfilecoll.GetList())

    # Process
    histos = drawer_book(options)
    drawer_project(tchain, histos, options)
    drawer_draw(histos, options)
    drawer_sitrep(histos, options)

# ______________________________________________________________________________
if __name__ == '__main__':

    # Setup argument parser
    parser = MyParser()
    options = parser.parse_args()

    # Call the main function
    main(options)
