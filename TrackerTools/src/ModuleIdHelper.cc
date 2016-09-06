#include "L1TTSimulations/TrackerTools/interface/ModuleIdHelper.h"

#include <cassert>

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"


uint32_t ModuleIdHelper::getModuleId(const TrackerTopology * tTopo, const DetId& detid) {
    uint32_t layer  = 999999;
    uint32_t ladder = 999999;
    uint32_t module = 999999;
    uint32_t type   = 999999;
    uint32_t ret    = 999999;

    if (detid.det() == DetId::Tracker) {
        //if (id.subdetId() == (int) PixelSubdetector::PixelBarrel) {
        //    PXBDetId pxbId(id);
        //    layer  = pxbId.layer();
        //    ladder = pxbId.ladder();
        //    module = pxbId.module();
        //
        //} else if (id.subdetId() == (int) PixelSubdetector::PixelEndcap) {
        //    PXFDetId pxfId(id);
        //    layer  = (pxfId.side() == 2) ? pxfId.disk() : pxfId.disk()+7;
        //    //ladder = pxfId.ring();
        //    ladder = pxfId.blade();
        //    module = pxfId.module();
        //}

        // Copied from Seb Viret
        //   https://github.com/sviret/HL_LHC/blob/62X_stuff/Extractors/RecoExtractor/src/StubExtractor.cc
        if ( detid.subdetId()==StripSubdetector::TOB )
        {
            layer  = static_cast<int>(tTopo->layer(detid))+4;
            ladder = static_cast<int>(tTopo->tobRod(detid));
            module = static_cast<int>(tTopo->module(detid));
            type   = static_cast<int>(tTopo->tobSide(detid)); // Tilt-/Tilt+/Flat <-> 1/2/3
        }
        else if ( detid.subdetId()==StripSubdetector::TID )
        {
            layer  = 10+static_cast<int>(tTopo->tidWheel(detid))+abs(2-static_cast<int>(tTopo->side(detid)))*7;
            ladder = static_cast<int>(tTopo->tidRing(detid));
            module = static_cast<int>(tTopo->module(detid));
            type   = 0;
        }

    }

    if (layer != 999999 && ladder != 999999 && module != 999999 && type != 999999) {
        ret = 10000*layer + 100*(ladder-1) + (module-1)/2;
    }

    return ret;
}
