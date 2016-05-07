#!/usr/bin/env python

import json
from math import pi, sqrt, atan2

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
ttmap = json.load(open("../data/trigger_sector_map.json"), object_pairs_hook=convert_key_to_int)
vertexmap = json.load(open("../data/module_vertices.json"), object_pairs_hook=convert_key_to_int)

# Find boundaries
def get_boundaries(tt, tt_moduleIds):
    phimin, phimax, etamin, etamax = get_parameter_space(tt)
    phimean = average([phimin, phimax])

    # Store min, max values for each layer
    phimins, phimaxs = [], []
    zmins, zmaxs = [], []  # z for barrel, r for endcap
    for i in xrange(23):
        phimins.append(phimean); phimaxs.append(phimean)
        zmins.append(9999.); zmaxs.append(-9999.)

    # Check vertex positions of each module
    for moduleId in tt_moduleIds:
        xyz = vertexmap[moduleId]
        assert(len(xyz) == 12)
        cyl = convert_xyz_to_cyl(xyz)

        lay = decodeLayer(moduleId)
        if 5 <= lay <= 10:
            # Barrel rectangle modules
            assert(xyz[0+0] == xyz[9+0])
            assert(xyz[0+1] == xyz[9+1])
            assert(xyz[0+2] == xyz[3+2])
            assert(xyz[3+0] == xyz[6+0])
            assert(xyz[3+1] == xyz[6+1])
            assert(xyz[6+2] == xyz[9+2])

            phi0, phi1 = phimean, phimean
            z0, z1 = 9999, -9999
            for i in xrange(4):
                phi = round(cyl[3*i+1],6)
                z = round(cyl[3*i+2],6)
                if deltaPhi(phi0,phimean) > deltaPhi(phi,phimean):
                    phi0 = phi
                if deltaPhi(phi1,phimean) < deltaPhi(phi,phimean):
                    phi1 = phi
                if z0 > z:
                    z0 = z
                if z1 < z:
                    z1 = z

            if deltaPhi(phimins[lay],phimean) > deltaPhi(phi0,phimean):
                phimins[lay] = phi0
            if deltaPhi(phimaxs[lay],phimean) < deltaPhi(phi1,phimean):
                phimaxs[lay] = phi1
            if zmins[lay] > z0:
                zmins[lay] = z0
            if zmaxs[lay] < z1:
                zmaxs[lay] = z1

        elif 11 <= lay <= 15 or 18 <= lay <= 22:
            # Endcap rectangle modules
            eps = 1.1e-4
            assert(abs((xyz[0+0]-xyz[3+0]) - (xyz[9+0]-xyz[6+0])) < eps)
            assert(abs((xyz[0+0]-xyz[9+0]) - (xyz[3+0]-xyz[6+0])) < eps)
            assert(abs((xyz[0+1]-xyz[3+1]) - (xyz[9+1]-xyz[6+1])) < eps)
            assert(abs((xyz[0+1]-xyz[9+1]) - (xyz[3+1]-xyz[6+1])) < eps)
            assert(xyz[0+2] == xyz[3+2])
            assert(xyz[6+2] == xyz[9+2])
            assert(xyz[0+2] == xyz[9+2])

            phi0, phi1 = phimean, phimean
            r0, r1 = 9999, -9999
            for i in xrange(4):
                phi = round(cyl[3*i+1],6)
                r = round(cyl[3*i+0],6)
                if deltaPhi(phi0,phimean) > deltaPhi(phi,phimean):
                    phi0 = phi
                if deltaPhi(phi1,phimean) < deltaPhi(phi,phimean):
                    phi1 = phi
                if r0 > r:
                    r0 = r
                if r1 < r:
                    r1 = r

            if deltaPhi(phimins[lay],phimean) > deltaPhi(phi0,phimean):
                phimins[lay] = phi0
            if deltaPhi(phimaxs[lay],phimean) < deltaPhi(phi1,phimean):
                phimaxs[lay] = phi1
            if zmins[lay] > r0:
                zmins[lay] = r0
            if zmaxs[lay] < r1:
                zmaxs[lay] = r1

        else:
            raise Exception("Unexpected moduleId: %i" % moduleId)

    return (phimins, phimaxs, zmins, zmaxs)


# ______________________________________________________________________________
# Main
writeme = []
for tt in xrange(48):
    if tt not in ttmap:
        raise Exception("Failed to retrieve from ttmap: %i" % tt)

    # Get the boundaries
    (phimins, phimaxs, zmins, zmaxs) = get_boundaries(tt, ttmap[tt])

    for i in xrange(23):
        if zmins[i] < 9999. or zmaxs[i] > -9999.:
            s = "{0},{1},{2:.6f},{3:.6f},{4:.4f},{5:.4f}".format(tt,i,phimins[i],phimaxs[i],zmins[i],zmaxs[i])
            #print s
            writeme.append(s)


# ______________________________________________________________________________
# Write .csv file
with open("../data/trigger_sector_boundaries.csv", "w") as f:
    writeme = ["tt/I,layer/I,phimin/D,phimax/D,zmin_cm/D,zmax_cm/D"] + writeme
    f.write("\n".join(writeme))

# Write .json file
mymap = {}
with open("../data/trigger_sector_boundaries.csv", "r") as f:
    for line in f:
        if not line[0].isdigit():
            continue
        values = line.split(",")
        assert(len(values) == 6)

        # Convert to int or float
        values = [float(x) if "." in x else int(x) for x in values]

        key = values[0]*100 + values[1]
        values = values[2:]
        mymap[key] = values

json.dump(mymap, open("../data/trigger_sector_boundaries.json", "w"), sort_keys=True)
