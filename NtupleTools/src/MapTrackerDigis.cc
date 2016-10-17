#include "L1TTSimulations/NtupleTools/interface/MapTrackerDigis.h"

void MapTrackerDigis::setup(const edm::Handle<collection_t>& handle) {
    if (handle.isValid()) {
        unsigned n = 0;
        for (digi_dsv_t::const_iterator itv = handle->begin(); itv != handle->end(); ++itv) {
            const DetId geoId(itv->detId());

            for (digi_ds_t::const_iterator it = itv->begin(); it != itv->end(); ++it) {
                //const Ref_Phase2TrackerDigi_ ref = edm::makeRefTo(handle, geoId, it);
                //mapping.insert(std::make_pair(ref, n));

                //const unsigned channel = it->channel();  // do not use as it contains also the overThreshold bit
                const unsigned channel = Phase2TrackerDigi::pixelToChannel(it->row(), it->column());
                const reference_t ref(geoId.rawId(), channel);
                const product_t prod = n;
                mapping.insert(std::make_pair(ref, prod));

                n++;
            }
        }
    }
}

int MapTrackerDigis::get(const reference_t& ref) const {
    return mapping.at(ref);
}
