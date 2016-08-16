#!/usr/bin/env python

from ROOT import TFile, TTree, TH1F, TH2F, TProfile, TF1, gROOT, gStyle, gPad
from math import sqrt, pi, sinh, atan2, tan, exp, log, degrees, radians
from itertools import izip, count
import numpy as np
from numpy import linalg as LA
from scipy import optimize
#import matplotlib.pyplot as plt
import random
random.seed(2016)

from incrementalstats import *

"""
# TODO
1. Derive pre-estimate
2. Allow missing layers
"""


# ______________________________________________________________________________
# Configurations
fname = "/cms/data/store/user/jiafulow/L1TrackTrigger/6_2_0_SLHC25p3/Demo_Emulator/stubs_tt27_300M_emu.root"
#fname = "root://cmsxrootd-site.fnal.gov//store/user/l1upgrades/SLHC/GEN/Demo_Emulator/stubs_tt27_300M_emu.root"
nentries = 100000
#nentries = 10000000
use_3D = 0
verbose = 1
make_plots = 1

nparameters2D = 2
nvariables2D = 6
nparameters3D = nparameters2D*2
nvariables3D = nvariables2D*2
cache_size = 10000

minInvPt = -1.0/3
maxInvPt = +1.0/3
minEta = 0.0
maxEta = 2.2/3

r_center = np.asarray([22.5913, 35.4772, 50.5402, 68.3101, 88.5002, 107.71])
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

    if use_3D:
        nvariables = nvariables3D
        nparameters = nparameters3D
    else:
        nvariables = nvariables2D
        nparameters = nparameters2D

    if verbose > 0:
        print "# Step 1"
        print "fname     : ", fname
        print "nentries  : ", nentries
        print "use_3D    : ", use_3D
        print "r_center  : ", r_center
        print "phi_center: ", phi_center
        print "eta_center: ", eta_center
        print

    # __________________________________________________________________________
    # Cache 10000 events

    stat_var1 = IncrementalStats(d=nvariables2D, cache_size=cache_size)
    stat_var2 = IncrementalStats(d=nvariables2D, cache_size=cache_size)

    for ievt, evt in enumerate(ttree):
        if ievt == nentries:
            break
        if stat_var1.count() == cache_size:
            break

        # Must have 1 particle
        if len(evt.genParts_pt) != 1:
            continue

        # Must have 6 stubs
        if len(evt.TTStubs_phi) != 6:
            continue

        # Sanity checks
        assert(len(evt.genParts_charge) == 1)
        assert(len(evt.genParts_pt) == 1)
        assert(len(evt.genParts_eta) == 1)
        assert(len(evt.TTStubs_phi) == 6)
        assert(len(evt.TTStubs_r) == 6)
        assert(len(evt.TTStubs_z) == 6)

        # Get track variables
        simInvPt = float(evt.genParts_charge[0])/evt.genParts_pt[0]
        simPhi = evt.genParts_phi[0]
        simEta = evt.genParts_eta[0]
        simCotTheta = sinh(simEta)
        simVz = evt.genParts_vz[0]

        try:
            theta = atan2(1.0, simCotTheta)
            eta = -log(tan(theta/2.0))
        except:
            print "Unexpected math error:", simEta, simCotTheta, simTanTheta, theta
            raise

        # Must satisfy invPt and eta ranges
        if not ((minInvPt <= simInvPt < maxInvPt) and (minEta <= simEta < maxEta)):
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

        # Apply R^3 correction
        r_over_two_rho_term = np.power(np.array(evt.TTStubs_r) * simC,3)
        r_over_two_rho_term *= 1.0/6.0
        r_over_two_rho_term_z = np.copy(r_over_two_rho_term)
        r_over_two_rho_term_z *= 1.0/simC
        variables1 -= r_over_two_rho_term
        variables2 -= r_over_two_rho_term_z

        stat_var1.add(variables1)
        stat_var2.add(variables2)
        continue

    if verbose > 0:
        print "count: ", stat_var1.count()
        print "mean : ", stat_var1.mean()
        print "var  : ", stat_var1.variance()
        print "cov  : ", stat_var1.covariance()
        print "q01  : ", stat_var1.quantile(p=0.01)
        print "q05  : ", stat_var1.quantile(p=0.05)
        print "q50  : ", stat_var1.quantile(p=0.50)
        print "q95  : ", stat_var1.quantile(p=0.95)
        print "q99  : ", stat_var1.quantile(p=0.99)
        print
        print "count: ", stat_var2.count()
        print "mean : ", stat_var2.mean()
        print "var  : ", stat_var2.variance()
        print "cov  : ", stat_var2.covariance()
        print "q01  : ", stat_var2.quantile(p=0.01)
        print "q05  : ", stat_var2.quantile(p=0.05)
        print "q50  : ", stat_var2.quantile(p=0.50)
        print "q95  : ", stat_var2.quantile(p=0.95)
        print "q99  : ", stat_var2.quantile(p=0.99)
        print

    # __________________________________________________________________________
    # Accumulate data

    #left_cuts_var1, right_cuts_var1 = stat_var1.quantile(p=0.05), stat_var1.quantile(p=0.95)
    #left_cuts_var2, right_cuts_var2 = stat_var2.quantile(p=0.05), stat_var2.quantile(p=0.95)
    left_cuts_var1, right_cuts_var1 = stat_var1.quantile(p=0.01), stat_var1.quantile(p=0.99)
    left_cuts_var2, right_cuts_var2 = stat_var2.quantile(p=0.01), stat_var2.quantile(p=0.99)

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
        if len(evt.TTStubs_phi) != 6:
            continue

        # Get track variables
        simInvPt = float(evt.genParts_charge[0])/evt.genParts_pt[0]
        simPhi = evt.genParts_phi[0]
        simEta = evt.genParts_eta[0]
        simCotTheta = sinh(simEta)
        simVz = evt.genParts_vz[0]
        parameters1 = np.array([simInvPt, simPhi])
        parameters2 = np.array([simCotTheta, simVz])

        # Must satisfy invPt and eta ranges
        if not ((minInvPt <= simInvPt < maxInvPt) and (minEta <= simEta < maxEta)):
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

        # Apply R^3 correction
        r_over_two_rho_term = np.power(np.array(evt.TTStubs_r) * simC,3)
        r_over_two_rho_term *= 1.0/6.0
        r_over_two_rho_term_z = np.copy(r_over_two_rho_term)
        r_over_two_rho_term_z *= 1.0/simC
        variables1 -= r_over_two_rho_term
        variables2 -= r_over_two_rho_term_z

        # Must satisfy variable ranges
        if not np.all([np.less_equal(left_cuts_var1, variables1), np.less(variables1, right_cuts_var1), np.less_equal(left_cuts_var2, variables2), np.less(variables2, right_cuts_var2)]):
            continue

        # Concatenate
        if use_3D:
            variables1 = np.concatenate((variables1, variables2))
            variables2 = variables1
            parameters1 = np.concatenate((parameters1, parameters2))
            parameters2 = parameters1

        stat_var1.add(variables1)
        stat_var2.add(variables2)
        stat_par1.add(parameters1)
        stat_par2.add(parameters2)

        if use_3D:
            data_var1.append(variables1)
            data_var2 = data_var1
            data_par1.append(parameters1)
            data_par2 = data_par1
        else:
            data_var1.append(variables1)
            data_var2.append(variables2)
            data_par1.append(parameters1)
            data_par2.append(parameters2)
        continue

    # Convert to numpy arrays
    #data_var1 = np.asarray(data_var1)
    #data_var2 = np.asarray(data_var2)
    #data_par1 = np.asarray(data_par1)
    #data_par2 = np.asarray(data_par2)

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

    if use_3D:
        nvariables = nvariables3D
        nparameters = nparameters3D
    else:
        nvariables = nvariables2D
        nparameters = nparameters2D

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
        principals1 = np.dot(V_var1, variables1)
        principals2 = np.dot(V_var2, variables2)

        stat_pc1.add(principals1)
        stat_pc2.add(principals2)
        data_pc1.append(principals1)
        data_pc2.append(principals2)

        # For D1 & D2
        stat_D1.add(principals1, covariables=parameters1)
        stat_D2.add(principals2, covariables=parameters2)
        continue

    # Solve for least squares D1 & D2
    #q_pc1, r_pc1 = LA.qr(stat_pc1.covariance())
    #d_pc1 = np.dot(LA.inv(r_pc1), np.dot(q_pc1.transpose(), stat_pc1.covariance()))
    #d_pc1 = LA.solve(stat_pc1.covariance(), stat_D1.covariance())
    soln_pc1 = LA.lstsq(stat_pc1.covariance(), stat_D1.covariance())
    soln_pc2 = LA.lstsq(stat_pc2.covariance(), stat_D2.covariance())
    d_pc1 = soln_pc1[0]
    d_pc2 = soln_pc2[0]
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
        #print "mean : ", stat_D1.mean()
        #print "var  : ", stat_D1.variance()
        print "cov  : ", stat_D1.covariance()
        print "D    : ", D_pc1
        print "DV   : ", D_var1
        print
        print "count: ", stat_D2.count()
        #print "mean : ", stat_D2.mean()
        #print "var  : ", stat_D2.variance()
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
def process_step2a(ntrials=2, do_trim=True):
    global stat_var1
    global stat_var2
    global stat_par1
    global stat_par2

    global w_var1
    global w_var2
    global V_var1
    global V_var2
    global D_var1
    global D_var2

    if use_3D:
        nvariables = nvariables3D
        nparameters = nparameters3D
    else:
        nvariables = nvariables2D
        nparameters = nparameters2D

    def weight_func(u):
        return None

    good_events = [True for i in xrange(len(data_var1))]

    for itrial in xrange(ntrials):
        if verbose > 0:
            print "# Step 2a"
            print "itrial: %i" % itrial

        stat_err1 = IncrementalStats(d=nparameters, cache_size=cache_size)
        stat_err2 = IncrementalStats(d=nparameters, cache_size=cache_size)

        for ievt, variables1, variables2, parameters1, parameters2 in izip(count(), data_var1, data_var2, data_par1, data_par2):
            if ievt == nentries:
                break
            if stat_err1.count() == cache_size:
                break

            # Skip events
            if not good_events[ievt]:
                continue

            # Fit parameters
            parameters_fit1 = np.dot(D_var1, variables1)
            parameters_fit2 = np.dot(D_var2, variables2)
            parameters_err1 = parameters_fit1 - parameters1
            parameters_err2 = parameters_fit2 - parameters2

            stat_err1.add(parameters_err1)
            stat_err2.add(parameters_err2)
            continue

        q25_err1 = stat_err1.quantile(p=0.25)
        q25_err2 = stat_err2.quantile(p=0.25)
        q75_err1 = stat_err1.quantile(p=0.75)
        q75_err2 = stat_err2.quantile(p=0.75)
        iqr_err1 = (q75_err1 - q25_err1)
        iqr_err2 = (q75_err2 - q25_err2)
        left_cuts_err1 = q25_err1 - iqr_err1 * 4
        left_cuts_err2 = q25_err2 - iqr_err2 * 4
        right_cuts_err1 = q75_err1 + iqr_err1 * 4
        right_cuts_err2 = q75_err2 + iqr_err2 * 4

        if verbose > 0:
            print "count : ", stat_err1.count()
            print "mean  : ", stat_err1.mean()
            print "std   : ", stat_err1.stdev()
            print "min   : ", stat_err1.minimum()
            print "max   : ", stat_err1.maximum()
            print "lcuts : ", left_cuts_err1
            print "rcuts : ", right_cuts_err1
            print
            print "count : ", stat_err2.count()
            print "mean  : ", stat_err2.mean()
            print "std   : ", stat_err2.stdev()
            print "min   : ", stat_err2.minimum()
            print "max   : ", stat_err2.maximum()
            print "lcuts : ", left_cuts_err2
            print "rcuts : ", right_cuts_err2
            print

        stat_var1 = IncrementalStats(d=nvariables)
        stat_var2 = IncrementalStats(d=nvariables)
        stat_par1 = IncrementalStats(d=nparameters)
        stat_par2 = IncrementalStats(d=nparameters)
        stat_pc1 = IncrementalStats(d=nvariables)
        stat_pc2 = IncrementalStats(d=nvariables)
        stat_D1 = IncrementalStats(d=nvariables, cov_d=nparameters)
        stat_D2 = IncrementalStats(d=nvariables, cov_d=nparameters)

        for ievt, variables1, variables2, parameters1, parameters2 in izip(count(), data_var1, data_var2, data_par1, data_par2):

            # Fit parameters
            parameters_fit1 = np.dot(D_var1, variables1)
            parameters_fit2 = np.dot(D_var2, variables2)
            parameters_err1 = parameters_fit1 - parameters1
            parameters_err2 = parameters_fit2 - parameters2
            w_err1 = weight_func(parameters_err1)
            w_err2 = weight_func(parameters_err2)

            # Simple trimming
            if do_trim:
                if not np.all([np.less_equal(left_cuts_err1, parameters_err1), np.less(parameters_err1, right_cuts_err1), np.less_equal(left_cuts_err2, parameters_err2), np.less(parameters_err2, right_cuts_err2)]):
                    good_events[ievt] = False

            # Skip events
            if not good_events[ievt]:
                continue

            stat_var1.add(variables1, weight=w_err1)
            stat_var2.add(variables2, weight=w_err2)
            stat_par1.add(parameters1, weight=w_err1)
            stat_par2.add(parameters2, weight=w_err2)

            # Principal components
            principals1 = np.dot(V_var1, variables1)
            principals2 = np.dot(V_var2, variables2)
            stat_pc1.add(principals1, weight=w_err1)
            stat_pc2.add(principals2, weight=w_err2)

            # For D1 & D2
            stat_D1.add(principals1, covariables=parameters1, weight=w_err1)
            stat_D2.add(principals2, covariables=parameters2, weight=w_err2)
            continue

        # Find eigenvectors
        w_var1, v_var1 = LA.eigh(stat_var1.covariance())
        w_var2, v_var2 = LA.eigh(stat_var2.covariance())
        V_var1 = v_var1.transpose()
        V_var2 = v_var2.transpose()

        # Solve for least squares D1 & D2
        soln_pc1 = LA.lstsq(stat_pc1.covariance(), stat_D1.covariance())
        soln_pc2 = LA.lstsq(stat_pc2.covariance(), stat_D2.covariance())
        d_pc1 = soln_pc1[0]
        d_pc2 = soln_pc2[0]
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
            #print "mean : ", stat_D1.mean()
            #print "var  : ", stat_D1.variance()
            print "cov  : ", stat_D1.covariance()
            print "D    : ", D_pc1
            print "DV   : ", D_var1
            print
            print "count: ", stat_D2.count()
            #print "mean : ", stat_D2.mean()
            #print "var  : ", stat_D2.variance()
            print "cov  : ", stat_D2.covariance()
            print "D    : ", D_pc2
            print "DV   : ", D_var2
            print
        continue
    return

# ______________________________________________________________________________
def process_step3():
    if use_3D:
        nvariables = nvariables3D
        nparameters = nparameters3D
    else:
        nvariables = nvariables2D
        nparameters = nparameters2D

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

    data_errPt = []

    for ievt, variables1, variables2, parameters1, parameters2 in izip(count(), data_var1, data_var2, data_par1, data_par2):
        # Principal components
        principals1 = np.dot(V_var1, variables1)
        principals2 = np.dot(V_var2, variables2)

        stat_pc1.add(principals1)
        stat_pc2.add(principals2)
        data_pc1.append(principals1)
        data_pc2.append(principals2)

        # Fit parameters
        parameters_fit1 = np.dot(D_var1, variables1)
        parameters_fit2 = np.dot(D_var2, variables2)
        parameters_err1 = parameters_fit1 - parameters1
        parameters_err2 = parameters_fit2 - parameters2

        stat_err1.add(parameters_err1)
        stat_err2.add(parameters_err2)
        data_err1.append(parameters_err1)
        data_err2.append(parameters_err2)

        errPt = (1.0/abs(parameters_fit1[0]) - 1.0/abs(parameters1[0])) * abs(parameters1[0])
        errInvPt = (parameters1[0] - parameters_fit1[0]) / abs(parameters1[0])
        data_errPt.append([errPt, errInvPt])
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

        myptbins = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 100.0, 125.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 2000.0, 7000.0]

        # Book histograms
        for i in xrange(nvariables3D):
            hname = "npc%i" % (i)
            htitle = "principal component %i" % (i)
            nbinsx, xmin, xmax = 1000, -7., 7.
            histos[hname] = TH1F(hname, ";"+htitle, nbinsx, xmin, xmax)

            htitle2, nbinsy, ymin, ymax = htitle, nbinsx, xmin, xmax
            for j in xrange(nparameters3D):
                hname = "npc%i_vs_par%i" % (i,j)
                if j == 0:
                    htitle = "q/p_{T} [1/GeV]"
                    nbinsx, xmin, xmax = 20, -0.5, 0.5
                elif j == 1:
                    htitle = "#phi [rad]"
                    nbinsx, xmin, xmax = 20, pi/4, pi/2
                elif j == 2:
                    htitle = "cot #theta"
                    nbinsx, xmin, xmax = 20, 0., 0.8
                elif j == 3:
                    htitle = "z_{0} [cm]"
                    nbinsx, xmin, xmax = 20, -15., 15.
                else:
                    htitle = ""
                    nbinsx, xmin, xmax = 20, -1., 1.
                histos[hname] = TH2F(hname, ";"+htitle+";"+htitle2, nbinsx, xmin, xmax, nbinsy, ymin, ymax)

        for i in xrange(nparameters3D + 2):
            hname = "err%i" % (i)
            if i == 0:
                htitle = "#Delta q/p_{T} [1/GeV]"
                nbinsx, xmin, xmax = 1000, -0.015, 0.015
            elif i == 1:
                htitle = "#Delta #phi [rad]"
                nbinsx, xmin, xmax = 1000, -0.005, 0.005
            elif i == 2:
                htitle = "#Delta cot #theta"
                nbinsx, xmin, xmax = 1000, -0.03, 0.03
            elif i == 3:
                htitle = "#Delta z_{0} [cm]"
                nbinsx, xmin, xmax = 1000, -1.2, 1.2
            elif i == 4:
                htitle = "#Delta (p_{T})/p_{T}"
                nbinsx, xmin, xmax = 1000, -0.2, 0.2
            elif i == 5:
                htitle = "#Delta (q/p_{T})*p_{T}"
                nbinsx, xmin, xmax = 1000, -0.2, 0.2
            else:
                htitle = ""
                nbinsx, xmin, xmax = 1000, -1., 1.
            histos[hname] = TH1F(hname, ";"+htitle, nbinsx, xmin, xmax)

            htitle2, nbinsy, ymin, ymax = htitle, nbinsx, xmin, xmax
            for j in xrange(nparameters3D):
                hname = "err%i_vs_par%i" % (i,j)
                if j == 0:
                    htitle = "q/p_{T} [1/GeV]"
                    nbinsx, xmin, xmax = 20, -0.5, 0.5
                elif j == 1:
                    htitle = "#phi [rad]"
                    nbinsx, xmin, xmax = 20, pi/4, pi/2
                elif j == 2:
                    htitle = "cot #theta"
                    nbinsx, xmin, xmax = 20, 0, 0.8
                elif j == 3:
                    htitle = "z_{0} [cm]"
                    nbinsx, xmin, xmax = 20, -15, 15
                else:
                    htitle = ""
                    nbinsx, xmin, xmax = 20, -1, 1
                histos[hname] = TH2F(hname, ";"+htitle+";"+htitle2, nbinsx, xmin, xmax, nbinsy, ymin, ymax)

            hname = "err%i_vs_pt" % (i)
            htitle = "p_{T} [GeV]"
            nbinsx, xmin, xmax = 20, 0, 200
            #histos[hname] = TH2F(hname, ";"+htitle+";"+htitle2, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
            histos[hname] = TH2F(hname, ";"+htitle+";"+htitle2, len(myptbins)-1, np.array(myptbins), nbinsy, ymin, ymax)
            hname = "err%i_vs_eta" % (i)
            htitle = "|#eta|"
            nbinsx, xmin, xmax = 25, 0, 2.5
            histos[hname] = TH2F(hname, ";"+htitle+";"+htitle2, nbinsx, xmin, xmax, nbinsy, ymin, ymax)


        # Fill histograms
        assert(len(data_par2) == len(data_par1))
        assert(len(data_pc1) == len(data_par1))
        assert(len(data_err1) == len(data_par1))

        for x1, x2, p1, p2 in izip(data_pc1, data_pc2, data_par1, data_par2):
            # Concatenate
            if not use_3D:
                x1 = [(x1[i] - stat_pc1.mean()[i])/sqrt(w_var1[i]) if abs(w_var1[i]) > 1e-14 else 0. for i in xrange(len(x1))]
                x2 = [(x2[i] - stat_pc2.mean()[i])/sqrt(w_var2[i]) if abs(w_var2[i]) > 1e-14 else 0. for i in xrange(len(x2))]
                x = np.concatenate((x1,x2))
                p = np.concatenate((p1,p2))
            else:
                x1 = [(x1[i] - stat_pc1.mean()[i])/sqrt(w_var1[i]) if abs(w_var1[i]) > 1e-14 else 0. for i in xrange(len(x1))]
                x = np.concatenate((x1,))
                p = np.concatenate((p1,))

            for i in xrange(nvariables3D):
                hname = "npc%i" % (i)
                histos[hname].Fill(x[i])

                for j in xrange(nparameters3D):
                    hname = "npc%i_vs_par%i" % (i,j)
                    histos[hname].Fill(p[j], x[i])

        for x1, x2, xErrPt, p1, p2 in izip(data_err1, data_err2, data_errPt, data_par1, data_par2):
            # Concatenate
            if not use_3D:
                x = np.concatenate((x1,x2,xErrPt))
                p = np.concatenate((p1,p2))
            else:
                x = np.concatenate((x1,xErrPt))
                p = np.concatenate((p1,))

            for i in xrange(nparameters3D + 2):
                hname = "err%i" % (i)
                histos[hname].Fill(x[i])

                for j in xrange(nparameters3D):
                    hname = "err%i_vs_par%i" % (i,j)
                    histos[hname].Fill(p[j], x[i])

                pt = 1.0/abs(p[0])
                theta = atan2(1.0, p[2])
                eta = -log(tan(theta/2.0))
                eta = abs(eta)
                hname = "err%i_vs_pt" % (i)
                histos[hname].Fill(pt, x[i])
                hname = "err%i_vs_eta" % (i)
                histos[hname].Fill(eta, x[i])

        # Write histograms
        outfile = TFile.Open("histos.root", "RECREATE")
        for hname, h in sorted(histos.iteritems()):
            h.Write()
        outfile.Close()

        # Print statistics
        printme = []
        for hname, h in sorted(histos.iteritems()):
            if h.ClassName() == "TH1F":
                h.Draw("hist")
                if h.Integral() > 0:
                    h.Fit("gaus","q")
                    h.fit = h.GetFunction("gaus")
                else:
                    h.fit = TF1("fa1", "gaus(0)")
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
    process_step2a()

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
