#!/usr/bin/env python

from rootdrawing import *
from helper import *
from math import pi, cos, sin, tan, sinh, asin, atan2, exp

# ______________________________________________________________________________
# Functions
quadsum = lambda x,y: sqrt(x*x + y*y)

# ______________________________________________________________________________
# Globals

# Load trigger tower map
ttmap = json.load(open("../data/trigger_sector_map_new.json"), object_pairs_hook=convert_key_to_int)

# Load module vertices
vertexmap = json.load(open("../data/module_vertices.json"), object_pairs_hook=convert_key_to_int)

moduleIds_set = set()
tpolylines_xy = {}
tpolylines_rz = {}
for moduleId, xyz in vertexmap.iteritems():
    if moduleId < 0:  continue
    moduleIds_set.add(moduleId)

    xy, rz = [], []
    for i in xrange(4):
        c = (xyz[3*i+0], xyz[3*i+1])
        if c not in xy:  xy.append(c)
        c = (xyz[3*i+2], quadsum(xyz[3*i], xyz[3*i+1]))
        if c not in rz:  rz.append(c)

    for coords, tpolylines in [(xy, tpolylines_xy), (rz, tpolylines_rz)]:
        x, y = [], []
        for c in coords:
            x.append(c[0])
            y.append(c[1])
        if len(xy) == 4:  # close the poly line in case it's made of 4 distinct points
            x.append(x[0])
            y.append(y[0])
        tpolylines[moduleId] = TPolyLine(len(x), array('d', x), array('d', y))

#count = 0
#for tt in xrange(48):
#    tt_moduleIds = ttmap[tt]
#    count += len(tt_moduleIds)
#
#print "Number of unique modules: %i" % len(moduleIds_set)
#print "Number of modules incl. duplicates: %i" % count

nbinsx, xmin, xmax = 100, -120., 120.
nbinsy, ymin, ymax = 100, -120., 120.
nbinsr, rmin, rmax = 100,    0., 120.
nbinsz, zmin, zmax = 100, -300., 300.

tline2 = TLine()
tline2.SetLineColor(kGray)
tline2.SetLineStyle(7)
tline2.SetLineWidth(1)

# ______________________________________________________________________________
# Drawer
def drawer_book(options):
    histos = {}

    hname = "xy1"
    histos[hname] = TH2F(hname, "; x [cm]; y [cm]", nbinsx, xmin, xmax, nbinsy, ymin, ymax)
    hname = "xy2"
    histos[hname] = TH2F(hname, "; x [cm]; y [cm]", nbinsx, xmin, xmax, nbinsy, ymin, ymax)
    hname = "rz"
    histos[hname] = TH2F(hname, "; z [cm]; r [cm]", nbinsz, zmin, zmax, nbinsy, rmin, rmax)

    # Style
    for hname, h in histos.iteritems():
        h.style = 0; h.logx = options.logx; h.logy = options.logy
    donotdelete.append(histos)
    return histos

# ______________________________________________________________________________
def drawer_project(tree, histos, options):
    return

# ______________________________________________________________________________
def drawer_draw(histos, options, debug=False):
    min_pt, max_vz, max_rho = options.min_pt, options.max_vz, options.max_rho
    tt = options.tt

    phimin, phimax, etamin, etamax = get_parameter_space(tt)
    cotmin = sinh(etamin)
    cotmax = sinh(etamax)

    mPtFactor = 0.3*3.8*1e-2/2.0
    invPt = 1.0/float(min_pt)

    # Magic numbers from tklayout
    #tklayout_phi_magic = pi/16
    #tklayout_z_magic = max_vz/max_rho

    # Magic numbers translated into rstar
    #rstar = sin(tklayout_phi_magic) / mPtFactor / invPt
    #rstar_z = max_vz/tklayout_z_magic
    #rstar = 63.4
    #rstar_z = max_rho
    rstar = 77.15
    rstar_z = 53.0

    if debug:  print "min_pt={0} max_vz={1} max_rho={2}".format(min_pt, max_vz, max_rho)
    if debug:  print "phimin={0:.4f} phimax={1:.4f} etamin={2:.4f} etamax={3:.4f} cotmin={4:.4f} cotmax={5:.4f}".format(phimin, phimax, etamin, etamax, cotmin, cotmax)
    if debug:  print "rstar={0:.4f} rstar_z={1:.4f}".format(rstar, rstar_z)

    # __________________________________________________________________________
    def traj_phimin_1(r):
        dphi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimin + dphi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_phimin_2(r):
        dphi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimin - dphi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_phimax_1(r):
        dphi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimax + dphi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_phimax_2(r):
        dphi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimax - dphi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_zmin_1(r):
        deltaZ = r * (cotmin + max_vz/rstar_z)
        z = -max_vz + deltaZ
        return z

    def traj_zmin_2(r):
        deltaZ = r * (cotmin - max_vz/rstar_z)
        z = +max_vz + deltaZ
        return z

    def traj_zmax_1(r):
        deltaZ = r * (cotmax + max_vz/rstar_z)
        z = -max_vz + deltaZ
        return z

    def traj_zmax_2(r):
        deltaZ = r * (cotmax - max_vz/rstar_z)
        z = +max_vz + deltaZ
        return z

    def traj_rmin_1(z):
        #if cotmax == 0:
        #    return 0
        r = (z + max_vz) / (cotmax + max_vz/rstar_z)
        return r

    def traj_rmin_2(z):
        #if cotmax == 0:
        #    return 0
        r = (z - max_vz) / (cotmax - max_vz/rstar_z)
        return r

    def traj_rmax_1(z):
        #if cotmin == 0:
        #    return 0
        r = (z + max_vz) / (cotmin + max_vz/rstar_z)
        return r

    def traj_rmax_2(z):
        #if cotmin == 0:
        #    return 0
        r = (z - max_vz) / (cotmin - max_vz/rstar_z)
        return r

    # __________________________________________________________________________
    def draw_phi_lines():
        for iphi in xrange(8):
            phi = -pi/2 + (2*pi/8) * iphi
            hyp = quadsum(xmax, ymax)
            #tline2.DrawLine(0, 0, hyp*cos(phi), hyp*sin(phi))
            if abs(hyp*cos(phi)) > xmax+0.1:
                hyp = xmax
            elif abs(hyp*sin(phi)) > ymax+0.1:
                hyp = quadsum(ymax/tan(phi), ymax)
            tline2.DrawLine(0, 0, hyp*cos(phi), hyp*sin(phi))

    def draw_eta_lines():
        for ieta in xrange(6+1):
            eta = -2.2 + (4.4/6) * ieta
            theta = 2.0 * atan2(exp(-eta),1)
            hyp = quadsum(zmax, rmax)
            #tline2.DrawLine(0, 0, hyp*cos(theta), hyp*sin(theta))
            if abs(hyp*cos(theta)) > zmax+0.1:
                hyp = zmax
            elif abs(hyp*sin(theta)) > rmax+0.1:
                hyp = quadsum(rmax/tan(theta), rmax)
            tline2.DrawLine(0, 0, hyp*cos(theta), hyp*sin(theta))

    # Make frames
    h = histos["xy1"]
    h.c1 = TCanvas("c_"+h.GetName(), "", 600, 600)
    h.SetStats(0); h.Draw()
    draw_phi_lines()

    h = histos["xy2"]
    h.c1 = TCanvas("c_"+h.GetName(), "", 600, 600)
    h.SetStats(0); h.Draw()
    draw_phi_lines()

    h = histos["rz"]
    h.c1 = TCanvas("c_"+h.GetName(), "", 800, 400)
    h.c1.SetLeftMargin(0.1); h.GetYaxis().SetTitleOffset(0.6)
    h.SetStats(0); h.Draw()
    draw_eta_lines()

    def style_tpolyline(l, module, colorized):
        l.SetLineStyle(1)
        l.SetLineWidth(2)
        if colorized:
            #if isPSModule(module):
            #    l.SetLineColor(palette[0])
            #else:
            #    l.SetLineColor(palette[1])
            l.SetLineColor(palette[3])
        else:
            l.SetLineColor(kGray)

    def doit():
        if 5 <= lay <= 10:  # Barrel
            histos["xy1"].c1.cd()
            tpolyline_xy.Draw()
            histos["rz"].c1.cd()
            tpolyline_rz.Draw()
        elif 11 <= lay <= 15 or 18 <= lay <= 22:  # Endcap
            histos["xy2"].c1.cd()
            tpolyline_xy.Draw()
            histos["rz"].c1.cd()
            tpolyline_rz.Draw()
        else:
            raise Exception("Unexpected moduleId: %i" % moduleId)

    # Draw polylines
    tpolylines = tpolylines_xy
    for moduleId in moduleIds_set:
        lay = decodeLayer(moduleId)
        tpolyline_xy = tpolylines_xy[moduleId]
        tpolyline_rz = tpolylines_rz[moduleId]
        style_tpolyline(tpolyline_xy, moduleId, 0)
        style_tpolyline(tpolyline_rz, moduleId, 0)
        doit()

    # Draw polylines for specific tower
    for moduleId in ttmap[tt]:
        lay = decodeLayer(moduleId)
        tpolyline_xy = tpolylines_xy[moduleId]
        tpolyline_rz = tpolylines_rz[moduleId]
        style_tpolyline(tpolyline_xy, moduleId, 1)
        style_tpolyline(tpolyline_rz, moduleId, 1)
        doit()

    # Make curves for 2 GeV muons
    def make_xy_graph(traj):
        n = 100
        xx, yy = [], []
        for i in xrange(n):
            r = rmin + (rmax - rmin) / float(n) * (i+0.5)
            phi = traj(r)
            xx.append(r * cos(phi))
            yy.append(r * sin(phi))
        gr = TGraph(n, array('d', xx), array('d', yy))
        return gr

    def make_rz_graph(traj):
        n = 100
        rr, zz = [], []
        for i in xrange(n):
            r = rmin + (rmax - rmin) / float(n) * (i+0.5)
            z = traj(r)
            rr.append(r)
            zz.append(z)
        gr = TGraph(n, array('d', zz), array('d', rr))
        return gr

    def make_rz_graph_using_z(traj):
        n = 100
        rr, zz = [], []
        for i in xrange(n):
            z = zmin + (zmax - zmin) / float(n) * (i+0.5)
            r = traj(z)
            rr.append(r)
            zz.append(z)
        gr = TGraph(n, array('d', zz), array('d', rr))
        return gr

    gr_phimin_1 = make_xy_graph(traj_phimin_1)
    gr_phimin_2 = make_xy_graph(traj_phimin_2)
    gr_phimax_1 = make_xy_graph(traj_phimax_1)
    gr_phimax_2 = make_xy_graph(traj_phimax_2)
    tgraphs = [gr_phimin_1, gr_phimin_2, gr_phimax_1, gr_phimax_2]

    gr_zmin_1 = make_rz_graph(traj_zmin_1)
    gr_zmin_2 = make_rz_graph(traj_zmin_2)
    gr_zmax_1 = make_rz_graph(traj_zmax_1)
    gr_zmax_2 = make_rz_graph(traj_zmax_2)
    #gr_zmin_1 = make_rz_graph_using_z(traj_rmin_1)
    #gr_zmin_2 = make_rz_graph_using_z(traj_rmin_2)
    #gr_zmax_1 = make_rz_graph_using_z(traj_rmax_1)
    #gr_zmax_2 = make_rz_graph_using_z(traj_rmax_2)
    tgraphs += [gr_zmin_1, gr_zmin_2, gr_zmax_1, gr_zmax_2]

    def style_tgraph(g, c):
        g.SetLineStyle(1)
        g.SetLineWidth(2)
        if c == 1:
            g.SetLineColor(palette[1])
        elif c == 2:
            g.SetLineColor(palette[0])

    for gr in [gr_phimin_1, gr_phimax_1, gr_zmin_1, gr_zmax_1]:
        style_tgraph(gr, 1)
    for gr in [gr_phimin_2, gr_phimax_2, gr_zmin_2, gr_zmax_2]:
        style_tgraph(gr, 2)

    histos["xy1"].c1.cd()
    for gr in [gr_phimin_1, gr_phimin_2, gr_phimax_1, gr_phimax_2]:
        gr.Draw("l")
    save(options.outdir, "xy1_tt%i" % options.tt, dot_pdf=False)

    histos["xy2"].c1.cd()
    for gr in [gr_phimin_1, gr_phimin_2, gr_phimax_1, gr_phimax_2]:
        gr.Draw("l")
    save(options.outdir, "xy2_tt%i" % options.tt, dot_pdf=False)

    histos["rz"].c1.cd()
    for gr in [gr_zmin_1, gr_zmin_2, gr_zmax_1, gr_zmax_2]:
        gr.Draw("l")
    save(options.outdir, "rz_tt%i" % options.tt, dot_pdf=False)

    donotdelete.append(tgraphs)
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
    #tchain.AddFileInfoList(options.tfilecoll.GetList())

    # Process
    histos = drawer_book(options)
    drawer_project(tchain, histos, options)
    drawer_draw(histos, options)
    drawer_sitrep(histos, options)

# ______________________________________________________________________________
if __name__ == '__main__':

    # Setup argument parser
    parser = MyParser()
    parser.add_argument("--min-pt", type=float, default=2., help="minimum track pT (default: %(default)s)")
    parser.add_argument("--max-vz", type=float, default=7., help="maximum vertex spread (default: %(default)s)")
    parser.add_argument("--max-rho", type=float, default=110., help="maximum tracker radius (default: %(default)s)")
    parser.add_argument("--tt", type=int, default=27, help="trigger tower (default: %(default)s)")
    options = parser.parse_args()

    # Call the main function
    main(options)

