#ifndef NtupleTools_MapTrackingParticles_h_
#define NtupleTools_MapTrackingParticles_h_

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include <map>


class MapTrackingParticles {
  public:
    typedef std::pair<unsigned, EncodedEventId> reference_t;
    typedef TrackingParticleRef                 product_t;
    typedef TrackingParticleCollection          collection_t;

    void setup(const edm::Handle<collection_t>& handle);

    unsigned size() const { return mapping.size(); }

    int get(const reference_t& ref) const;

  private:
    std::map<reference_t, product_t> mapping;
};

#endif
