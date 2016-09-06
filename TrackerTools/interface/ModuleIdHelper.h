#ifndef TrackerTools_ModuleIdHelper_h_
#define TrackerTools_ModuleIdHelper_h_

#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
//#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
//#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"


class TrackerTopology;

class ModuleIdHelper {
  public:
    static uint32_t getModuleId(const TrackerTopology * tTopo, const DetId& detid);
};

#endif
