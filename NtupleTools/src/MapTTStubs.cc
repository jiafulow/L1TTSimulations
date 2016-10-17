#include "L1TTSimulations/NtupleTools/interface/MapTTStubs.h"

void MapTTStubs::setup(const edm::Handle<collection_t>& handle) {
    if (handle.isValid()) {
        unsigned n = 0;
        for (ttstub_dsv_t::const_iterator itv = handle->begin(); itv != handle->end(); ++itv) {
            for (ttstub_ds_t::const_iterator it = itv->begin(); it != itv->end(); ++it) {
                const reference_t ref = edmNew::makeRefTo(handle, it);
                const product_t prod = n;
                mapping.insert(std::make_pair(ref, prod));

                n++;
            }
        }
    }
}

int MapTTStubs::get(const reference_t& ref) const {
    return mapping.at(ref);
}
