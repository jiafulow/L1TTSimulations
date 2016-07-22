#!/usr/bin/env python

from ROOT import TFile, TTree, TH1F, gROOT, gStyle, gPad
from math import sqrt, pi, sinh, atan2
from itertools import izip, count
import numpy as np
from numpy import linalg as LA
from scipy import optimize
#import matplotlib.pyplot as plt
import random

from incrementalstats import *

"""
# TODO
1. Derive pre-estimate
2. Plot err as a function of parameters
3. Allow missing layers
"""


# ______________________________________________________________________________
# Configurations
fname = "/cms/data/store/user/jiafulow/L1TrackTrigger/6_2_0_SLHC25p3/Demo_Emulator/stubs_tt27_300M_emu.root"
#fname = "root://cmsxrootd-site.fnal.gov//store/user/l1upgrades/SLHC/GEN/Demo_Emulator/stubs_tt27_300M_emu.root"
nentries = 100000
#nentries = 10000000
verbose = 1
make_plots = 1

nparameters = 2  # 2D view
nvariables = 6   # 2D view
cache_size = 10000

minInvPt = -1.0/3
maxInvPt = +1.0/3

r_center = np.array([22.5913, 35.4772, 50.5402, 68.3101, 88.5002, 107.71])
phi_center = pi*3/8
eta_center = 2.2/6
cot_center = sinh(eta_center)

# Data
stat_var1 = None
stat_var2 = None
stat_par1 = None
stat_par2 = None
data_var1 = None
data_var2 = None
data_par1 = None
data_par2 = None

# Results
w_var1 = None
w_var2 = None
V_var1 = None
V_var2 = None
D_var1 = None
D_var2 = None


# ______________________________________________________________________________
def process_step1():
    global stat_var1
    global stat_var2
    global stat_par1
    global stat_par2
    global data_var1
    global data_var2
    global data_par1
    global data_par2

    if verbose > 0:
        print "# Step 1"
        print "r_center  : ", r_center
        print "phi_center: ", phi_center
        print "eta_center: ", eta_center
        print

    # __________________________________________________________________________
    # Cache 10000 events

    stat_var1 = IncrementalStats(d=nvariables, cache_size=cache_size)
    stat_var2 = IncrementalStats(d=nvariables, cache_size=cache_size)

    for ievt, evt in enumerate(ttree):
        if ievt == min(cache_size, nentries):
            break

        # Must have 1 particle
        if len(evt.genParts_pt) != 1:
            continue

        # Must have 6 stubs
        if len(evt.TTStubs_phi) != nvariables:
            continue

        # Sanity checks
        assert(len(evt.genParts_charge) == 1)
        assert(len(evt.genParts_pt) == 1)
        assert(len(evt.genParts_eta) == 1)
        assert(len(evt.TTStubs_phi) == nvariables)
        assert(len(evt.TTStubs_r) == nvariables)
        assert(len(evt.TTStubs_z) == nvariables)

        # Get track variables
        simInvPt = float(evt.genParts_charge[0])/evt.genParts_pt[0]
        simCotTheta = sinh(evt.genParts_eta[0])

        # Must satisfy invPt range
        if not (minInvPt <= simInvPt < maxInvPt):
            continue

        # Get stub variables
        variables1 = np.array(evt.TTStubs_phi)
        variables2 = np.array(evt.TTStubs_z)
        variables3 = np.array(evt.TTStubs_r) - r_center  # deltaR

        # Apply z smearing
        variables2[0] += 0.5*random.uniform(-0.15, 0.15)
        variables2[1] += 0.5*random.uniform(-0.15, 0.15)
        variables2[2] += 0.5*random.uniform(-0.15, 0.15)
        variables2[3] += 0.5*random.uniform(-5, 5)
        variables2[4] += 0.5*random.uniform(-5, 5)
        variables2[5] += 0.5*random.uniform(-5, 5)

        # Apply deltaR correction
        simC = -0.5 * (0.003 * 3.8 * simInvPt)  # 1/(2 x radius of curvature)
        simT = simCotTheta
        variables1 -= simC * variables3
        variables2 -= simT * variables3

        stat_var1.add(variables1)
        stat_var2.add(variables2)
        continue

    if verbose > 0:
        print "count: ", stat_var1.count()
        print "mean : ", stat_var1.mean()
        print "var  : ", stat_var1.variance()
        print "cov  : ", stat_var1.covariance()
        print "q05  : ", stat_var1.quantile(p=0.05)
        print "q50  : ", stat_var1.quantile(p=0.50)
        print "q95  : ", stat_var1.quantile(p=0.95)
        print
        print "count: ", stat_var2.count()
        print "mean : ", stat_var2.mean()
        print "var  : ", stat_var2.variance()
        print "cov  : ", stat_var2.covariance()
        print "q05  : ", stat_var2.quantile(p=0.05)
        print "q50  : ", stat_var2.quantile(p=0.50)
        print "q95  : ", stat_var2.quantile(p=0.95)
        print

    # __________________________________________________________________________
    # Accumulate data

    min_var1, max_var1 = stat_var1.quantile(p=0.05), stat_var1.quantile(p=0.95)
    min_var2, max_var2 = stat_var2.quantile(p=0.05), stat_var2.quantile(p=0.95)

    stat_var1 = IncrementalStats(d=nvariables)
    stat_var2 = IncrementalStats(d=nvariables)
    stat_par1 = IncrementalStats(d=nparameters)
    stat_par2 = IncrementalStats(d=nparameters)

    data_var1 = []
    data_var2 = []
    data_par1 = []
    data_par2 = []

    for ievt, evt in enumerate(ttree):
        if ievt == nentries:
            break

        # Must have 1 particle
        if len(evt.genParts_pt) != 1:
            continue

        # Must have 6 stubs
        if len(evt.TTStubs_phi) != nvariables:
            continue

        # Get track variables
        simInvPt = float(evt.genParts_charge[0])/evt.genParts_pt[0]
        simPhi = evt.genParts_phi[0]
        simCotTheta = sinh(evt.genParts_eta[0])
        simVz = evt.genParts_vz[0]
        parameters1 = np.array([simInvPt, simPhi])
        parameters2 = np.array([simCotTheta, simVz])

        # Must satisfy invPt range
        if not (minInvPt <= simInvPt < maxInvPt):
            continue

        # Get stub variables
        variables1 = np.array(evt.TTStubs_phi)
        variables2 = np.array(evt.TTStubs_z)
        variables3 = np.array(evt.TTStubs_r) - r_center  # deltaR

        # Apply z smearing
        variables2[0] += 0.5*random.uniform(-0.15, 0.15)
        variables2[1] += 0.5*random.uniform(-0.15, 0.15)
        variables2[2] += 0.5*random.uniform(-0.15, 0.15)
        variables2[3] += 0.5*random.uniform(-5, 5)
        variables2[4] += 0.5*random.uniform(-5, 5)
        variables2[5] += 0.5*random.uniform(-5, 5)

        # Apply deltaR correction
        simC = -0.5 * (0.003 * 3.8 * simInvPt)  # 1/(2 x radius of curvature)
        simT = simCotTheta
        variables1 -= simC * variables3
        variables2 -= simT * variables3

        # Must satisfy variable ranges
        if not np.all([np.less_equal(min_var1, variables1), np.less(variables1, max_var1), np.less_equal(min_var2, variables2), np.less(variables2, max_var2)]):
            continue

        stat_var1.add(variables1)
        stat_var2.add(variables2)
        stat_par1.add(parameters1)
        stat_par2.add(parameters2)

        data_var1.append(variables1)
        data_var2.append(variables2)
        data_par1.append(parameters1)
        data_par2.append(parameters2)
        continue

    if verbose > 0:
        print "count: ", stat_var1.count()
        print "mean : ", stat_var1.mean()
        print "var  : ", stat_var1.variance()
        print "cov  : ", stat_var1.covariance()
        print
        print "count: ", stat_var2.count()
        print "mean : ", stat_var2.mean()
        print "var  : ", stat_var2.variance()
        print "cov  : ", stat_var2.covariance()
        print
        print "count: ", stat_par1.count()
        print "mean : ", stat_par1.mean()
        print "var  : ", stat_par1.variance()
        print "cov  : ", stat_par1.covariance()
        print
        print "count: ", stat_par2.count()
        print "mean : ", stat_par2.mean()
        print "var  : ", stat_par2.variance()
        print "cov  : ", stat_par2.covariance()
        print
    return

# ______________________________________________________________________________
def process_step2():
    global w_var1
    global w_var2
    global V_var1
    global V_var2
    global D_var1
    global D_var2

    # Find eigenvectors
    w_var1, v_var1 = LA.eigh(stat_var1.covariance())
    w_var2, v_var2 = LA.eigh(stat_var2.covariance())
    V_var1 = v_var1.transpose()
    V_var2 = v_var2.transpose()

    if verbose > 0:
        print "# Step 2"
        print "w: ", w_var1
        print "v: ", v_var1
        print
        print "w: ", w_var2
        print "v: ", v_var2
        print
        # Verify
        #print "verify: ", np.dot(stat_var1.covariance(), v_var1[:, 0]) - w_var1[0] * v_var1[:, 0]
        #print "verify: ", np.dot(stat_var2.covariance(), v_var2[:, 0]) - w_var2[0] * v_var2[:, 0]
        #print

    # __________________________________________________________________________
    # Find solutions

    stat_pc1 = IncrementalStats(d=nvariables)
    stat_pc2 = IncrementalStats(d=nvariables)
    data_pc1 = []
    data_pc2 = []

    stat_D1 = IncrementalStats(d=nvariables, cov_d=nparameters)
    stat_D2 = IncrementalStats(d=nvariables, cov_d=nparameters)

    for ievt, variables1, variables2, parameters1, parameters2 in izip(count(), data_var1, data_var2, data_par1, data_par2):
        # Principal components
        principals1 = np.dot(V_var1, variables1 - stat_var1.mean())
        principals2 = np.dot(V_var2, variables2 - stat_var2.mean())

        stat_pc1.add(principals1)
        stat_pc2.add(principals2)
        data_pc1.append(principals1)
        data_pc2.append(principals2)

        # For D1 & D2
        stat_D1.add(principals1, covariables=parameters1)
        stat_D2.add(principals2, covariables=parameters2)
        continue

    #q_pc1, r_pc1 = LA.qr(stat_pc1.covariance())
    #d_pc1 = np.dot(LA.inv(r_pc1), np.dot(q_pc1.transpose(), stat_pc1.covariance()))
    d_pc1 = LA.solve(stat_pc1.covariance(), stat_D1.covariance())
    d_pc2 = LA.solve(stat_pc2.covariance(), stat_D2.covariance())
    D_pc1 = d_pc1.transpose()
    D_pc2 = d_pc2.transpose()
    D_var1 = np.dot(D_pc1, V_var1)
    D_var2 = np.dot(D_pc2, V_var2)

    if verbose > 0:
        print "count: ", stat_pc1.count()
        print "mean : ", stat_pc1.mean()
        print "var  : ", stat_pc1.variance()
        print "cov  : ", stat_pc1.covariance()
        print
        print "count: ", stat_pc2.count()
        print "mean : ", stat_pc2.mean()
        print "var  : ", stat_pc2.variance()
        print "cov  : ", stat_pc2.covariance()
        print
        print "count: ", stat_D1.count()
        print "mean : ", stat_D1.mean()
        print "var  : ", stat_D1.variance()
        print "cov  : ", stat_D1.covariance()
        print "D    : ", D_pc1
        print "DV   : ", D_var1
        print
        print "count: ", stat_D2.count()
        print "mean : ", stat_D2.mean()
        print "var  : ", stat_D2.variance()
        print "cov  : ", stat_D2.covariance()
        print "D    : ", D_pc2
        print "DV   : ", D_var2
        print
        # Verify
        #print "verify: ", np.dot(stat_pc1.covariance(), soln_D1)
        #print "verify: ", np.dot(stat_pc2.covariance(), soln_D2)
        #print
    return

# ______________________________________________________________________________
def process_step3():
    if verbose > 0:
        print "# Step 3"

    # __________________________________________________________________________
    # Evaluate fit errors

    stat_pc1 = IncrementalStats(d=nvariables, cache_size=cache_size)
    stat_pc2 = IncrementalStats(d=nvariables, cache_size=cache_size)
    data_pc1 = []
    data_pc2 = []

    stat_err1 = IncrementalStats(d=nparameters, cache_size=cache_size)
    stat_err2 = IncrementalStats(d=nparameters, cache_size=cache_size)
    data_err1 = []
    data_err2 = []

    for ievt, variables1, variables2, parameters1, parameters2 in izip(count(), data_var1, data_var2, data_par1, data_par2):
        # Principal components
        principals1 = np.dot(V_var1, variables1 - stat_var1.mean())
        principals2 = np.dot(V_var2, variables2 - stat_var2.mean())

        stat_pc1.add(principals1)
        stat_pc2.add(principals2)
        data_pc1.append(principals1)
        data_pc2.append(principals2)

        # Fit parameters
        #parameters_fit1 = np.dot(D_var1, variables1 - stat_var1.mean())
        #parameters_fit2 = np.dot(D_var2, variables2 - stat_var2.mean())
        parameters_fit1 = np.dot(D_var1, variables1)
        parameters_fit2 = np.dot(D_var2, variables2)
        parameters_err1 = parameters_fit1 - parameters1
        parameters_err2 = parameters_fit2 - parameters2

        stat_err1.add(parameters_err1)
        stat_err2.add(parameters_err2)
        data_err1.append(parameters_err1)
        data_err2.append(parameters_err2)
        continue

    if verbose > 0:
        print "count: ", stat_pc1.count()
        print "mean : ", stat_pc1.mean()
        print "var  : ", stat_pc1.variance()
        print "cov  : ", stat_pc1.covariance()
        print "q05  : ", stat_pc1.quantile(p=0.05)
        print "q50  : ", stat_pc1.quantile(p=0.50)
        print "q95  : ", stat_pc1.quantile(p=0.95)
        print
        print "count: ", stat_pc2.count()
        print "mean : ", stat_pc2.mean()
        print "var  : ", stat_pc2.variance()
        print "cov  : ", stat_pc2.covariance()
        print "q05  : ", stat_pc2.quantile(p=0.05)
        print "q50  : ", stat_pc2.quantile(p=0.50)
        print "q95  : ", stat_pc2.quantile(p=0.95)
        print
        print "count: ", stat_err1.count()
        print "mean : ", stat_err1.mean()
        print "var  : ", stat_err1.variance()
        print "cov  : ", stat_err1.covariance()
        print "q05  : ", stat_err1.quantile(p=0.05)
        print "q50  : ", stat_err1.quantile(p=0.50)
        print "q95  : ", stat_err1.quantile(p=0.95)
        print
        print "count: ", stat_err2.count()
        print "mean : ", stat_err2.mean()
        print "var  : ", stat_err2.variance()
        print "cov  : ", stat_err2.covariance()
        print "q05  : ", stat_err2.quantile(p=0.05)
        print "q50  : ", stat_err2.quantile(p=0.50)
        print "q95  : ", stat_err2.quantile(p=0.95)
        print

    # __________________________________________________________________________
    # Make plots
    if make_plots:
        histos = {}

        # Book histograms
        for i in xrange(nvariables*2):
            hname = "npc%i" % (i)
            nbinsx, xmin, xmax = 1000, -7., 7.
            histos[hname] = TH1F(hname, ";"+hname, nbinsx, xmin, xmax)

        for i in xrange(nparameters*2):
            hname = "err%i" % (i)
            if i == 0:
                nbinsx, xmin, xmax = 1000, -0.02, 0.02
            elif i == 1:
                nbinsx, xmin, xmax = 1000, -0.005, 0.005
            elif i == 2:
                nbinsx, xmin, xmax = 1000, -0.02, 0.02
            elif i == 3:
                nbinsx, xmin, xmax = 1000, -0.5, 0.5
            else:
                nbinsx, xmin, xmax = 1000, -2., 2.
            histos[hname] = TH1F(hname, ";"+hname, nbinsx, xmin, xmax)

        # Fill histograms
        for x in data_pc1:
            for i in xrange(nvariables):
                xx = x[i]/sqrt(w_var1[i]) if abs(w_var1[i]) > 1e-14 else 0.
                hname = "npc%i" % (i)
                histos[hname].Fill(xx)

        for x in data_pc2:
            for i in xrange(nvariables):
                xx = x[i]/sqrt(w_var2[i]) if abs(w_var2[i]) > 1e-14 else 0.
                hname = "npc%i" % (nvariables+i)
                histos[hname].Fill(xx)

        for x in data_err1:
            for i in xrange(nparameters):
                xx = x[i]
                hname = "err%i" % (i)
                histos[hname].Fill(xx)

        for x in data_err2:
            for i in xrange(nparameters):
                xx = x[i]
                hname = "err%i" % (nparameters+i)
                histos[hname].Fill(xx)

        # Write histograms
        printme = []
        for hname, h in sorted(histos.iteritems()):
            h.Draw("hist")
            h.Fit("gaus","q")
            h.fit = h.GetFunction("gaus")
            imgname = "%s.png" % hname
            gPad.Print(imgname)

            s = "%s, %f, %f, %f, %f" % (hname, h.GetMean(), h.GetRMS(), h.fit.GetParameter(1), h.fit.GetParameter(2))
            printme.append(s)
        print '\n'.join(printme)

    return

# ______________________________________________________________________________
def main():
    # 1st step: loop over events to get rough estimates of all the statistics,
    #   reject outliers, accumulate data
    process_step1()

    # 2nd step: compute the coefficients
    process_step2()

    # 3rd step: evaluate with the coefficients
    process_step3()

    return

# ______________________________________________________________________________
if __name__ == '__main__':

    tfile = TFile.Open(fname)
    ttree = tfile.Get("ntupler/tree")
    gROOT.LoadMacro("tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle()")
    gROOT.SetBatch(True)
    gStyle.SetOptStat(111110)

    main()
