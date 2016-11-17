#!/usr/bin/env python

from ROOT import TChain, TH1F, TColor, gPad, gStyle, gSystem, gROOT
from itertools import izip, count
from array import array
import argparse

# ______________________________________________________________________________
# Drawer
def drawer_book(histos, options):
    hname = "nroads_per_road"
    histos[hname] = TH1F(hname, "; # roads/road", 11, -0.5, 10.5)

    hname = "nroads_per_event"
    histos[hname] = TH1F(hname, "; # roads/tower/BX", 600, 0., 600.)

    hname = "ncombinations_per_road"
    histos[hname] = TH1F(hname, "; # combinations/road", 400, 0., 400.)

    hname = "ncombinations_per_event"
    histos[hname] = TH1F(hname, "; # combinations/tower/BX", 1600, 0., 1600.)

    hname = "nreadouts_per_road"
    histos[hname] = TH1F(hname, "; # readouts/road", 11, -0.5, 10.5)

    hname = "nreadouts_per_event"
    histos[hname] = TH1F(hname, "; # readouts/tower/BX", 600, 0., 600.)

    # Get colors
    if options.is_merged:
      col  = TColor.GetColor("#1A1AE3")  # merged
      fcol = TColor.GetColor("#9999FB")  # merged
    else:
      col  = TColor.GetColor("#e31a1c")  # unmerged
      fcol = TColor.GetColor("#fb9a99")  # unmerged

    # Style
    for hname, h in histos.iteritems():
        h.SetLineWidth(2); h.SetMarkerSize(0)
        h.SetLineColor(col); h.SetFillColor(fcol)
        if h.ClassName() == "TH1F":
            binwidth = (h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin())/h.GetNbinsX()
            h.SetYTitle("Entries / %.1f" % binwidth)
    return

# ______________________________________________________________________________
def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("AMTTRoads_patternRef"         , 1)
    tree.SetBranchStatus("AMTTRoads_nstubs"             , 1)
    tree.SetBranchStatus("AMTTRoads_stubRefs"           , 1)
    tree.SetBranchStatus("AMTTRoads_patternInvPt"       , 1)
    tree.SetBranchStatus("AMTTRoads_superstripIds"      , 1)
    tree.SetBranchStatus("AMTTRoads_superstripIdsUnited", 1)


    # __________________________________________________________________________
    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break
        if (ievt % 1000 == 0):  print "Processing event: %i" % ievt

        if not options.is_merged:
            evt.AMTTRoads_superstripIdsUnited = evt.AMTTRoads_superstripIds

        if options.verbose:
            print ".. %i nroads: %i" % (ievt, len(evt.AMTTRoads_superstripIdsUnited))

        nroads_per_event        = 0
        ncombinations_per_event = 0
        nreadouts_per_event        = 0

        # Loop over roads
        for iroad, patternRef, nstubs, stubRefs, patternInvPt, superstripIds in izip(
            count(), evt.AMTTRoads_patternRef, evt.AMTTRoads_nstubs, evt.AMTTRoads_stubRefs, evt.AMTTRoads_patternInvPt, evt.AMTTRoads_superstripIdsUnited
        ):

            if patternRef >= options.maxPatterns:
                continue

            if not options.is_merged:
                superstripIds = [[x] for x in superstripIds]

            if options.verbose:
                for ilayer, ilayer_stubRefs, ilayer_superstripIds in izip(count(), stubRefs, superstripIds):
                    print ".... %i nstubs: %i nsuperstrips: %i" % (ilayer, len(ilayer_stubRefs), len(ilayer_superstripIds))

            nroads_per_road = 1

            ncombinations_per_road = 1
            for ilayer, ilayer_stubRefs in izip(count(), stubRefs):
                nstubs_per_layer = len(ilayer_stubRefs)
                if nstubs_per_layer != 0:
                    ncombinations_per_road *= nstubs_per_layer

            nreadouts_per_road = 0
            for ilayer, ilayer_superstripIds in izip(count(), superstripIds):
                nsuperstrips_per_layer = len(ilayer_superstripIds)
                if nreadouts_per_road < nsuperstrips_per_layer:
                    nreadouts_per_road = nsuperstrips_per_layer

            # Fill histos
            histos["nroads_per_road"].Fill(nroads_per_road)
            histos["ncombinations_per_road"].Fill(ncombinations_per_road)
            histos["nreadouts_per_road"].Fill(nreadouts_per_road)

            nroads_per_event += 1
            ncombinations_per_event += ncombinations_per_road
            nreadouts_per_event += nreadouts_per_road
            continue

        # Fill histos
        histos["nroads_per_event"].Fill(nroads_per_event)
        histos["ncombinations_per_event"].Fill(ncombinations_per_event)
        histos["nreadouts_per_event"].Fill(nreadouts_per_event)
        continue

    tree.SetBranchStatus("*", 1)
    return

# ______________________________________________________________________________
def drawer_draw(histos, options):
    def displayQuantiles(h, in_quantiles=[0.95,0.99], scalebox=(1.,1.)):
        # Display one-sided confidence intervals, a.k.a quantiles
        n = len(in_quantiles)
        in_quantiles = array('d', in_quantiles)
        quantiles = array('d', [0. for i in xrange(n)])
        h.GetQuantiles(n, quantiles, in_quantiles)

        gPad.Modified(); gPad.Update()
        ps = h.FindObject("stats")
        ps.SetName("mystats")

        newX1NDC = ps.GetX2NDC() - (ps.GetX2NDC() - ps.GetX1NDC()) * scalebox[0]
        newY1NDC = ps.GetY2NDC() - ((ps.GetY2NDC() - ps.GetY1NDC()) / 5 * (5 + n)) * scalebox[1]
        ps.SetX1NDC(newX1NDC)
        ps.SetY1NDC(newY1NDC)

        for iq, q in enumerate(in_quantiles):
            ps.AddText("%i%% CI = %6.4g" % (int(q*100), quantiles[iq]))
        h.stats = [h.GetMean()] + quantiles.tolist()

        h.SetStats(0)
        #gPad.Modified(); gPad.Update()
        ps.Draw()

    for hname, h in histos.iteritems():
        if h.ClassName() == "TH1F":
            if options.logy:
                h.SetMaximum(h.GetMaximum() * 14); h.SetMinimum(0.5)
            else:
                h.SetMaximum(h.GetMaximum() * 1.4); h.SetMinimum(0.)
            h.SetStats(1); h.Draw("hist")
            gPad.SetLogy(options.logy)
            displayQuantiles(h)

            #CMS_options.label()
            #save(options.outdir, hname)
            gPad.Print("%s/%s_%s.png" % (options.outdir, options.label, hname))
    return

# ______________________________________________________________________________
def drawer_sitrep(histos, options):
    print "--- SITREP ---------------------------------------------------------"
    print
    return


# ______________________________________________________________________________
# Main function
def main(histos, options):
    tchain = TChain("ntupler/tree", "")
    tchain.AddFile(options.filename)

    drawer_book(histos, options)
    drawer_project(tchain, histos, options)
    drawer_draw(histos, options)
    drawer_sitrep(histos, options)

# ______________________________________________________________________________
if __name__ == '__main__':

    gROOT.SetBatch(True)

    gROOT.LoadMacro("tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle();")

    gStyle.SetOptStat(111110)
    gStyle.SetStatX(0.94)
    gStyle.SetStatY(0.93)
    gStyle.SetStatH(0.30)
    gStyle.SetStatW(0.28)

    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--program', type=int, default=1, help='an integer')
    options = parser.parse_args()
    options.nentries = 10000
    options.verbose = 1
    options.logy = True
    options.outdir = "figures_runRoadMerging/"

    # Create outdir if necessary
    if not options.outdir.endswith("/"):
        options.outdir += "/"
    if gSystem.AccessPathName(options.outdir):
        gSystem.mkdir(options.outdir)

    # Hardcoded configurations
    programs = {}

    # label, filename, is_merged, maxPatterns
    programs[0] = ("TTbar_PU140", "../dataFiles/roads_TTbar_PU140.root", False, 5532635)
    programs[1] = ("m8_TTbar_PU140", "../dataFiles/m8_roads_TTbar_PU140.root", True, 1129402)
    programs[2] = ("TTbar_PU200", "../dataFiles/roads_TTbar_PU200.root", False, 5532635)
    programs[3] = ("m8_TTbar_PU200", "../dataFiles/m8_roads_TTbar_PU200.root", True, 1129402)

    options.label, options.filename, options.is_merged, options.maxPatterns = programs[options.program]

    histos = {}

    main(histos, options)
