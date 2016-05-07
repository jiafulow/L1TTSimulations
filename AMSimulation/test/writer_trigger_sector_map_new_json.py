#!/usr/bin/env python

import json
from math import pi, sqrt, sin, sinh, asin, atan2

# tklayout version:
#   https://github.com/tkLayout/tkLayout/blob/master/src/AnalyzerVisitors/TriggerProcessorBandwidth.cpp

# ______________________________________________________________________________
# Functions
convert_key_to_int = lambda pairs: dict([(int(k),v) for (k,v) in pairs])

quadsum = lambda x,y: sqrt(x*x + y*y)

average = lambda x: sum(x, 0.0) / len(x)

def convert_xyz_to_cyl(xyz):
    cyl = []
    for i in xrange(4):
        r = quadsum(xyz[3*i], xyz[3*i+1])
        phi = atan2(xyz[3*i+1], xyz[3*i])
        z = xyz[3*i+2]
        cyl += [r, phi, z]
    return cyl

def deltaPhi(phi1, phi2):
    result = phi1 - phi2
    while result > pi:  result -= 2*pi
    while result <= -pi:  result += 2*pi
    return result

# Get layer number from moduleId
def decodeLayer(moduleId):
    return int(moduleId / 10000)

# Get the parameter space of a trigger tower
def get_parameter_space(tt):
    ieta = tt/8
    iphi = tt%8
    etamin = -2.2 + (4.4/6) * ieta
    etamax = -2.2 + (4.4/6) * (ieta+1)
    phimin = -pi/2 + (2*pi/8) * iphi
    phimax = -pi/2 + (2*pi/8) * (iphi+1)
    if iphi >= 6:
        phimin -= 2*pi
        phimax -= 2*pi
    return (phimin, phimax, etamin, etamax)


# ______________________________________________________________________________
# Load existing maps
#ttmap = json.load(open("../data/trigger_sector_map.json"), object_pairs_hook=convert_key_to_int)
vertexmap = json.load(open("../data/module_vertices.json"), object_pairs_hook=convert_key_to_int)

vertexmap_cyl = {}
moduleIds_set = set()
min_rho, max_rho = 9999, -9999
for moduleId, xyz in vertexmap.iteritems():
    cyl = convert_xyz_to_cyl(xyz)
    vertexmap_cyl[moduleId] = cyl
    if moduleId > 0:
        moduleIds_set.add(moduleId)
    if min_rho > average(cyl[0::3]):
        min_rho = average(cyl[0::3])
    if max_rho < average(cyl[0::3]):
        max_rho = average(cyl[0::3])
print "min_rho, max_rho = %.4f, %.4f" % (min_rho, max_rho)

# Get the trajectories in the physical space of a trigger tower
def get_trajectories(tt, min_pt=2., max_vz=7., max_rho=110., debug=False):
    phimin, phimax, etamin, etamax = get_parameter_space(tt)
    cotmin = sinh(etamin)
    cotmax = sinh(etamax)

    mPtFactor = 0.3*3.8*1e-2/2.0
    invPt = 1.0/float(min_pt)

    # Magic numbers from tklayout
    #tklayout_phi_magic = pi/24
    tklayout_phi_magic = pi/16
    tklayout_z_magic = max_vz/max_rho

    # Magic numbers translated into rstar
    rstar = sin(tklayout_phi_magic) / mPtFactor / invPt
    rstar_z = max_vz/tklayout_z_magic

    if debug:  print "min_pt={0} max_vz={1} max_rho={2}".format(min_pt, max_vz, max_rho)
    if debug:  print "phimin={0:.4f} phimax={1:.4f} etamin={2:.4f} etamax={3:.4f} cotmin={4:.4f} cotmax={5:.4f}".format(phimin, phimax, etamin, etamax, cotmin, cotmax)

    # Define the trajectories
    def traj_phimin(r):
        deltaPhi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimin + deltaPhi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_phimax(r):
        deltaPhi = - asin(mPtFactor * (r-rstar) * invPt)
        phi = phimax - deltaPhi
        if phi > +pi:  phi -= 2*pi
        if phi < -pi:  phi += 2*pi
        return phi

    def traj_zmin(r):
        #deltaZ = (1.0 / (mPtFactor * invPt) * asin(mPtFactor * r * invPt)) * cotmin
        #deltaZ = r * (cotmin + tklayout_z_magic)
        deltaZ = r * (cotmin + max_vz/rstar_z)
        return -max_vz + deltaZ

    def traj_zmax(r):
        #deltaZ = (1.0 / (mPtFactor * invPt) * asin(mPtFactor * r * invPt)) * cotmax
        #deltaZ = r * (cotmax - tklayout_z_magic)
        deltaZ = r * (cotmax - max_vz/rstar_z)
        return +max_vz + deltaZ

    def traj_rmin(z):
        if cotmax == 0:
            return 0
        #r = (z - max_vz) / (cotmax - tklayout_z_magic)
        r = (z - max_vz) / (cotmax - max_vz/rstar_z)
        return r

    def traj_rmax(z):
        if cotmin == 0:
            return 0
        #r = (z + max_vz) / (cotmin + tklayout_z_magic)
        r = (z + max_vz) / (cotmin + max_vz/rstar_z)
        return r

    return (traj_phimin, traj_phimax, traj_zmin, traj_zmax, traj_rmin, traj_rmax)

# Filter the modules based on the trajectories
def filter_modules(trajs, trajsL, trajsR, moduleIds_set, debug=False):
    results = []

    traj_phimin, traj_phimax, traj_zmin, traj_zmax, traj_rmin, traj_rmax = trajs
    #trajL_phimin, trajL_phimax, trajL_zmin, trajL_zmax, trajL_rmin, trajL_rmax = trajsL
    #trajR_phimin, trajR_phimax, trajR_zmin, trajR_zmax, trajR_rmin, trajR_rmax = trajsR

    def is_phi_good(phis, phimin, phimax):
        return any(deltaPhi(phimin, phimin) <= deltaPhi(phi, phimin) <= deltaPhi(phimax, phimin) for phi in phis)

    def is_z_good(zees, zmin, zmax):
        return any(zmin <= zee <= zmax for zee in zees)

    for moduleId in moduleIds_set:
        cyl = vertexmap_cyl[moduleId]

        lay = decodeLayer(moduleId)
        if 5 <= lay <= 10:  # Barrel
            r = average(cyl[0::3])
            phimin = traj_phimin(r)
            phimax = traj_phimax(r)
            #if (deltaPhi(phimin, trajL_phimax(r)) > 0):
            #    phimin = trajL_phimax(r)
            #if (deltaPhi(trajR_phimin(r), phimax) > 0):
            #    phimax = trajR_phimin(r)
            zmin = traj_zmin(r)
            zmax = traj_zmax(r)

            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f}".format(moduleId, r, phimin, phimax, zmin, zmax)
            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, r, *(cyl[1::3] + cyl[2::3]))

            assert(deltaPhi(phimax, phimin) > 0)
            assert(zmax - zmin > 0)

            if is_phi_good(cyl[1::3], phimin, phimax) and is_z_good(cyl[2::3], zmin, zmax):
                results.append(moduleId)

        elif 11 <= lay <= 15 or 18 <= lay <= 22:  # Endcap
            z = average(cyl[2::3])
            if z > 0:
                rmin = traj_rmin(z)
                rmax = traj_rmax(z)
            else:
                rmin = traj_rmax(z)
                rmax = traj_rmin(z)
            r = average(cyl[0::3])
            phimin = traj_phimin(r)
            phimax = traj_phimax(r)

            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f}".format(moduleId, z, phimin, phimax, rmin, rmax)
            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, z, *(cyl[1::3] + cyl[0::3]))

            assert(deltaPhi(phimax, phimin) > 0)
            assert((16 <= tt <= 31) or (rmax - rmin > 0))

            if is_phi_good(cyl[1::3], phimin, phimax) and is_z_good(cyl[0::3], rmin, rmax):
                results.append(moduleId)

    return set(results)


# ______________________________________________________________________________
# Main
writeme = []
count = 0
for tt in xrange(48):
    # Get the trajectories
    trajs = get_trajectories(tt)
    trajsL = get_trajectories((tt/8)*8+((tt-1)%8))
    trajsR = get_trajectories((tt/8)*8+((tt+1)%8))

    # Get the moduleIds
    tt_moduleIds = filter_modules(trajs, trajsL, trajsR, moduleIds_set)
    count += len(tt_moduleIds)

    ieta = tt/8
    iphi = tt%8
    s = ",".join("{0}".format(n) for n in [ieta+1,iphi+1]+sorted(tt_moduleIds))
    #print s
    writeme.append(s)

print "Number of unique modules: %i" % len(moduleIds_set)
print "Number of modules incl. duplicates: %i" % count


# ______________________________________________________________________________
# Write .csv file
with open("../data/trigger_sector_map_new.csv", "w") as f:
    writeme = ["eta_idx, phi_idx, module_list"] + writeme
    f.write("\n".join(writeme))

# Write .json file
mymap = {}
with open("../data/trigger_sector_map_new.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")

        # Convert to int
        values = [int(x) for x in values]

        key = (values[0]-1)*8 + (values[1]-1)
        values = sorted(values[2:])
        mymap[key] = values

assert(len(mymap) == 6*8)

json.dump(mymap, open("../data/trigger_sector_map_new.json", "w"), sort_keys=True)
