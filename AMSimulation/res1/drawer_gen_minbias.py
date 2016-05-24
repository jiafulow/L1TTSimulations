#!/usr/bin/env python

from rootdrawing import *
from math import sqrt, sinh

col  = TColor.GetColor("#e31a1c")  # nu140
fcol = TColor.GetColor("#fb9a99")  # nu140

pt_vec = [2,3,5,10]

# ______________________________________________________________________________
# Drawer
def drawer_book(options):
    histos = {}

    hname = "trkParts_pt"
    histos[hname] = TH1F(hname, "; trkPart p_{T} [GeV]", 200, 0, 100)

    hname = "trkParts_invPt"
    histos[hname] = TH1F(hname, "; trkPart q/p_{T} [1/GeV]", 200, -5, 5)

    for pt in pt_vec:
        hname = "trkParts_eta_pt%i" % pt
        histos[hname] = TH1F(hname, "; trkPart #eta", 200, -5, 5)

        hname = "trkParts_cottheta_pt%i" % pt
        histos[hname] = TH1F(hname, "; trkPart cot #theta", 200, -75, 75)

    hname = "trkParts_phi"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} #phi [rad]", 200, -3.2, 3.2)

    hname = "trkParts_vx"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} vertex x [cm]", 200, -0.01, 0.01)

    hname = "trkParts_vy"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} vertex y [cm]", 200, -0.01, 0.01)

    hname = "trkParts_vxy"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} vertex perp [cm]", 200, -0.01, 0.01)

    hname = "trkParts_vz"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} vertex z [cm]", 200, -30, 30)

    hname = "trkParts_charge"
    histos[hname] = TH1F(hname, "; trkPart_{pt>2} charge", 5, -2.5, 2.5)

    hname = "ntrkParts_pt0p5"
    histos[hname] = TH1F(hname, "; # trkParts {pt>0.5}", 250, 0, 2500)

    hname = "ntrkParts_pt1"
    histos[hname] = TH1F(hname, "; # trkParts {pt>1}", 250, 0, 1000)

    hname = "ntrkParts_pt2"
    histos[hname] = TH1F(hname, "; # trkParts {pt>2}", 250, 0, 250)

    hname = "ntrkParts_pt3"
    histos[hname] = TH1F(hname, "; # trkParts {pt>3}", 250, 0, 250)

    hname = "ntrkParts_pt5"
    histos[hname] = TH1F(hname, "; # trkParts {pt>5}", 250, 0, 250)

    hname = "ntrkParts_pt10"
    histos[hname] = TH1F(hname, "; # trkParts {pt>10}", 250, 0, 250)

    hname = "nvertices"
    histos[hname] = TH1F(hname, "; # vertices", 250, 0, 250)

    # Style
    for hname, h in histos.iteritems():
        h.style = 0; h.logx = options.logx; h.logy = options.logy
    donotdelete.append(histos)
    return histos

# ______________________________________________________________________________
def drawer_project(tree, histos, options):
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("trkParts_pt"     , 1)
    tree.SetBranchStatus("trkParts_eta"    , 1)
    tree.SetBranchStatus("trkParts_phi"    , 1)
    tree.SetBranchStatus("trkParts_vx"     , 1)
    tree.SetBranchStatus("trkParts_vy"     , 1)
    tree.SetBranchStatus("trkParts_vz"     , 1)
    tree.SetBranchStatus("trkParts_charge" , 1)
    tree.SetBranchStatus("trkParts_pdgId"  , 1)
    tree.SetBranchStatus("trkParts_primary", 1)
    tree.SetBranchStatus("gen_nPV"         , 1)

    # Loop over events
    for ievt, evt in enumerate(tree):
        if (ievt == options.nentries):  break
        if (ievt % 10 == 0):  print "Processing event: %i" % ievt

        nparts = evt.trkParts_pt.size()
        pt_map = {}

        for ipart in xrange(nparts):
            primary = evt.trkParts_primary[ipart]
            charge  = evt.trkParts_charge [ipart]

            if not (charge!=0 and primary):
                continue

            pt      = evt.trkParts_pt     [ipart]
            eta     = evt.trkParts_eta    [ipart]
            phi     = evt.trkParts_phi    [ipart]
            vx      = evt.trkParts_vx     [ipart]
            vy      = evt.trkParts_vy     [ipart]
            vz      = evt.trkParts_vz     [ipart]
            charge  = evt.trkParts_charge [ipart]
            pdgId   = evt.trkParts_pdgId  [ipart]
            cottheta = sinh(eta)

            if pt > 2:  # relax eta cut
                histos["trkParts_eta_pt2"].Fill(eta)
                histos["trkParts_cottheta_pt2"].Fill(cottheta)
            if pt > 3:  # relax eta cut
                histos["trkParts_eta_pt3"].Fill(eta)
                histos["trkParts_cottheta_pt3"].Fill(cottheta)
            if pt > 5:  # relax eta cut
                histos["trkParts_eta_pt5"].Fill(eta)
                histos["trkParts_cottheta_pt5"].Fill(cottheta)
            if pt > 10:  # relax eta cut
                histos["trkParts_eta_pt10"].Fill(eta)
                histos["trkParts_cottheta_pt10"].Fill(cottheta)

            if abs(eta) < 2.4:  # relax pt cut
                histos["trkParts_pt"     ].Fill(pt)
                histos["trkParts_invPt"  ].Fill(float(charge)/pt)

                pt_map[ipart] = pt

            if not (pt > 2 and abs(eta) < 2.4):
                histos["trkParts_phi"    ].Fill(phi)
                histos["trkParts_vx"     ].Fill(vx)
                histos["trkParts_vy"     ].Fill(vy)
                histos["trkParts_vxy"    ].Fill(sqrt(vx*vx+vy*vy))
                histos["trkParts_vz"     ].Fill(vz)
                histos["trkParts_charge" ].Fill(charge)

        histos["ntrkParts_pt0p5" ].Fill(count_if(pt_map.values(), lambda x: x>0.5))
        histos["ntrkParts_pt1"   ].Fill(count_if(pt_map.values(), lambda x: x>1))
        histos["ntrkParts_pt2"   ].Fill(count_if(pt_map.values(), lambda x: x>2))
        histos["ntrkParts_pt3"   ].Fill(count_if(pt_map.values(), lambda x: x>3))
        histos["ntrkParts_pt5"   ].Fill(count_if(pt_map.values(), lambda x: x>5))
        histos["ntrkParts_pt10"  ].Fill(count_if(pt_map.values(), lambda x: x>10))

        histos["nvertices"       ].Fill(evt.gen_nPV)

    tree.SetBranchStatus("*", 1)
    return

# ______________________________________________________________________________
def drawer_draw(histos, options):
    for hname, h in histos.iteritems():

        if h.ClassName() == "TH1F":
            if True:
                # Draw TH1 plots
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

    # Specialized
    if True:
        for hvar in ["eta", "cottheta"]:
            pt = pt_vec[0]
            hname = "trkParts_%s_pt%i" % (hvar, pt)
            h = histos[hname]
            h.SetMinimum(0)
            h.SetStats(0); h.Draw("hist")

            moveLegend(tlegend,0.20,0.78,0.55,0.92); tlegend.Clear()
            for i, pt in enumerate(pt_vec):
                hname = "trkParts_%s_pt%i" % (hvar, pt)
                h = histos[hname]
                h.SetLineColor(palette[i])
                h.Draw("hist same")
                tlegend.AddEntry(h, "p_{T} > %i GeV" % pt, "l")
            tlegend.Draw()

            CMS_label()
            save(options.outdir, "trkParts_%s_overlay" % (hvar))

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
