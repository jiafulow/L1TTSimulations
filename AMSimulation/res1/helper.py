from math import pi

# ______________________________________________________________________________
# Encoding

def halfStripRound(x):
    p = 10.
    return int(floor((x*2)*p + 0.5)/p)

def decodeLayer(moduleId):
    return int(moduleId / 10000)

def decodeLadder(moduleId):
    return int(moduleId / 100) % 100

def decodeModule(moduleId):
    return int(moduleId) % 100

def isPSModule(moduleId):
    lay = decodeLayer(moduleId)
    if 5 <= lay <= 7:
        return True
    lad = decodeLadder(moduleId)
    if 11 <= lay <= 22 and lad <= 8:
        return True
    return False

def isBarrelModule(moduleId):
    lay = decodeLayer(moduleId)
    if 5 <= lay <= 10:
        return True
    return False

def isEndcapModule(moduleId):
    lay = decodeLayer(moduleId)
    if 11 <= lay <= 15 or 18 <= lay <= 22:
        return True
    return False

# ______________________________________________________________________________
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

