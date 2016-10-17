#ifndef NtupleTools_MapTrackerDigis_h_
#define NtupleTools_MapTrackerDigis_h_

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include <map>


class MapTrackerDigis {
  public:
    typedef Phase2TrackerDigi                     digi_t;
    typedef edm::DetSetVector<digi_t>             digi_dsv_t;
    typedef edm::DetSet<digi_t>                   digi_ds_t;
    typedef Ref_Phase2TrackerDigi_                digi_ref_t;

    typedef std::pair<unsigned, unsigned>         reference_t;
    typedef unsigned                              product_t;
    typedef digi_dsv_t                            collection_t;

    void setup(const edm::Handle<collection_t>& handle);

    unsigned size() const { return mapping.size(); }

    int get(const reference_t& ref) const;

  private:
    std::map<reference_t, product_t> mapping;
};

#endif
