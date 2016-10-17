#include "L1TTSimulations/NtupleTools/interface/MapTrackingParticles.h"

void MapTrackingParticles::setup(const edm::Handle<collection_t>& handle) {
    if (handle.isValid()) {
        for (unsigned itp=0; itp<handle->size(); ++itp) {
            const product_t prod(handle, itp);
            for (TrackingParticle::g4t_iterator itrk = prod->g4Track_begin(); itrk != prod->g4Track_end(); ++itrk) {
                const reference_t ref(itrk->trackId(), prod->eventId());
                mapping.insert(std::make_pair(ref, prod));
            }
        }
    }
}

int MapTrackingParticles::get(const reference_t& ref) const {
    std::map<reference_t, product_t>::const_iterator found = mapping.find(ref);
    if (found != mapping.end()) {
        return found->second.key();
    }
    return -1;
}
