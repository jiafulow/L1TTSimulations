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

    // offsets_2023tilted[layer-5][type-1]
    static const int offsets_2023tilted[6][3] = {
      {0, 11+ 7, 11},
      {0, 12+11, 12},
      {0, 13+15, 13},
      {0, 0, 0},
      {0, 0, 0},
      {0, 0, 0},
    };

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

        if (detid.subdetId() == StripSubdetector::TOB) {  // Phase2 Outer Tracker Barrel
            layer  = 4+tTopo->layer(detid);   // increasing r
            ladder = tTopo->tobRod(detid);    // increasing abs(z) (rings) or phi (barrel)
            module = tTopo->module(detid);    // increasing phi (rings) or z (bareel)
            type   = tTopo->tobSide(detid);   // 1=rings- 2=rings+ 3=barrel0

            if (type == 1 || type == 2) {
              std::swap(ladder, module);
            }
            if (type == 1 || type == 2 || type == 3) {
              assert((layer-5) < 6 && (type-1) < 3);
              module += offsets_2023tilted[layer-5][type-1];
            }

        } else if (detid.subdetId() == StripSubdetector::TID) {  // Phase2 Outer Tracker Endcap
            layer  = 10+tTopo->layer(detid);  // increasing abs(z)
            ladder = tTopo->tidRing(detid);   // increasing r
            module = tTopo->module(detid);    // increasing phi
            type   = tTopo->tidSide(detid);   // 1=-ve 2=+ve

            if (type == 1) {
              layer += 7;
            }
        }
    }

    if (layer != 999999 && ladder != 999999 && module != 999999 && type != 999999) {
        ret = 10000*layer + 100*(ladder-1) + (module-1);
    }

    return ret;
}
