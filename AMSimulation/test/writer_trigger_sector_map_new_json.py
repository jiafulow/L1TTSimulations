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
min_rho_1, max_rho_1 = 9999., -9999.
for moduleId, xyz in vertexmap.iteritems():
    if moduleId < 0:  continue
    moduleIds_set.add(moduleId)

    cyl = convert_xyz_to_cyl(xyz)
    vertexmap_cyl[moduleId] = cyl

    if min_rho_1 > average(cyl[0::3]):
        min_rho_1 = average(cyl[0::3])
    if max_rho_1 < average(cyl[0::3]):
        max_rho_1 = average(cyl[0::3])
#print "min_rho_1, max_rho_1 = %.4f, %.4f" % (min_rho_1, max_rho_1)

# Get the trajectories in the physical space of a trigger tower
def get_trajectories(tt, min_pt=2., max_vz=7., max_rho=110., debug=False):
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

    # Define the trajectories
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

    return (traj_phimin_1, traj_phimin_2, traj_phimax_1, traj_phimax_2,
        traj_zmin_1, traj_zmin_2, traj_zmax_1, traj_zmax_2,
        traj_rmin_1, traj_rmin_2, traj_rmax_1, traj_rmax_2)

# Filter the modules based on the trajectories
def filter_modules(tt, trajs, moduleIds_set, debug=False):
    results = []

    (traj_phimin_1, traj_phimin_2, traj_phimax_1, traj_phimax_2,
        traj_zmin_1, traj_zmin_2, traj_zmax_1, traj_zmax_2,
        traj_rmin_1, traj_rmin_2, traj_rmax_1, traj_rmax_2) = trajs

    def is_phi_good(phis, phimin, phimax):
        return any(deltaPhi(phimin, phimin) <= deltaPhi(phi, phimin) <= deltaPhi(phimax, phimin) for phi in phis)

    def is_z_good(zees, zmin, zmax):
        return any(zmin <= zee <= zmax for zee in zees)

    for moduleId in moduleIds_set:
        cyl = vertexmap_cyl[moduleId]

        lay = decodeLayer(moduleId)
        if 5 <= lay <= 10:  # Barrel
            r = average(cyl[0::3])
            r_phimin_1 = traj_phimin_1(r)
            r_phimin_2 = traj_phimin_2(r)
            r_phimax_1 = traj_phimax_1(r)
            r_phimax_2 = traj_phimax_2(r)
            r_zmin_1 = traj_zmin_1(r)
            r_zmin_2 = traj_zmin_2(r)
            r_zmax_1 = traj_zmax_1(r)
            r_zmax_2 = traj_zmax_2(r)

            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, r, r_phimin_1, r_phimin_2, r_phimax_1, r_phimax_2, r_zmin_1, r_zmin_2, r_zmax_1, r_zmax_2)
            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, r, *(cyl[1::3] + cyl[2::3]))

            assert(deltaPhi(r_phimax_1, r_phimin_1) > 0)
            assert(deltaPhi(r_phimax_2, r_phimin_2) > 0)
            assert(r_zmax_1 - r_zmin_1 > 0)
            assert(r_zmax_2 - r_zmin_2 > 0)

            #if (is_phi_good(cyl[1::3], r_phimin_1, r_phimax_1) or is_phi_good(cyl[1::3], r_phimin_2, r_phimax_2)) and \
            #    (is_z_good(cyl[2::3], r_zmin_1, r_zmax_1) or is_z_good(cyl[2::3], r_zmin_2, r_zmax_2)):
            #    results.append(moduleId)

            # Decide to use whether 1 or 2
            if deltaPhi(r_phimax_1, r_phimin_1) >= deltaPhi(r_phimax_2, r_phimin_1):
                r_phimax = r_phimax_1
            else:
                r_phimax = r_phimax_2
            if deltaPhi(r_phimin_1, r_phimin_1) <= deltaPhi(r_phimin_2, r_phimin_1):
                r_phimin = r_phimin_1
            else:
                r_phimin = r_phimin_2
            r_zmax = max(r_zmax_1, r_zmax_2)
            r_zmin = min(r_zmin_1, r_zmin_2)

            if is_phi_good(cyl[1::3], r_phimin, r_phimax) and is_z_good(cyl[2::3], r_zmin, r_zmax):
                results.append(moduleId)

        elif 11 <= lay <= 15 or 18 <= lay <= 22:  # Endcap
            if 11 <= lay <= 15 and tt <  24:  continue  # ignore opposite endcap
            if 18 <= lay <= 22 and tt >= 24:  continue  # ignore opposite endcap

            z = average(cyl[2::3])
            if z > 0:
                z_rmin_1 = traj_rmin_1(z)
                z_rmin_2 = traj_rmin_2(z)
                z_rmax_1 = traj_rmax_1(z)
                z_rmax_2 = traj_rmax_2(z)
            else:
                z_rmin_1 = traj_rmax_1(z)
                z_rmin_2 = traj_rmax_2(z)
                z_rmax_1 = traj_rmin_1(z)
                z_rmax_2 = traj_rmin_2(z)
            r = average(cyl[0::3])
            r_phimin_1 = traj_phimin_1(r)
            r_phimin_2 = traj_phimin_2(r)
            r_phimax_1 = traj_phimax_1(r)
            r_phimax_2 = traj_phimax_2(r)

            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, z, r_phimin_1, r_phimin_2, r_phimax_1, r_phimax_2, z_rmin_1, z_rmin_2, z_rmax_1, z_rmax_2)
            if debug:  print "{0} {1:.4f} {2:.4f} {3:.4f} {4:.4f} {5:.4f} {6:.4f} {7:.4f} {8:.4f} {9:.4f}".format(moduleId, z, *(cyl[1::3] + cyl[0::3]))

            assert(deltaPhi(r_phimax_1, r_phimin_1) > 0)
            assert(deltaPhi(r_phimax_2, r_phimin_2) > 0)
            assert((z_rmax_1 < 0) or (z_rmax_1 - z_rmin_1 > 0))
            assert((z_rmax_2 < 0) or (z_rmax_2 - z_rmin_2 > 0))

            #if (is_phi_good(cyl[1::3], r_phimin_1, r_phimax_1) or is_phi_good(cyl[1::3], r_phimin_2, r_phimax_2)) and \
            #    (is_z_good(cyl[0::3], z_rmin_1, z_rmax_1) or is_z_good(cyl[0::3], z_rmin_2, z_rmax_2)):
            #    results.append(moduleId)

            # Decide to use whether 1 or 2
            if deltaPhi(r_phimax_1, r_phimin_1) >= deltaPhi(r_phimax_2, r_phimin_1):
                r_phimax = r_phimax_1
            else:
                r_phimax = r_phimax_2
            if deltaPhi(r_phimin_1, r_phimin_1) <= deltaPhi(r_phimin_2, r_phimin_1):
                r_phimin = r_phimin_1
            else:
                r_phimin = r_phimin_2
            safe_r = lambda r: r if r > 0 else 9999.
            z_rmax = max(safe_r(z_rmax_1), safe_r(z_rmax_2))
            z_rmin = min(safe_r(z_rmin_1), safe_r(z_rmin_2))

            if is_phi_good(cyl[1::3], r_phimin, r_phimax) and is_z_good(cyl[0::3], z_rmin, z_rmax):
                results.append(moduleId)

        else:
            raise Exception("Unexpected moduleId: %i" % moduleId)

    return set(results)


# ______________________________________________________________________________
# Main
writeme = []
count = 0
for tt in xrange(48):
    # Get the trajectories
    trajs = get_trajectories(tt)

    # Get the moduleIds
    tt_moduleIds = filter_modules(tt, trajs, moduleIds_set)
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
