#ifndef NtupleTools_MapTTStubs_h_
#define NtupleTools_MapTTStubs_h_

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include <map>


class MapTTStubs {
  public:
    typedef Ref_Phase2TrackerDigi_                digi_ref_t;
    typedef TTStub<digi_ref_t>                    ttstub_t;
    typedef edmNew::DetSetVector<ttstub_t>        ttstub_dsv_t;
    typedef edmNew::DetSet<ttstub_t>              ttstub_ds_t;
    typedef edm::Ref<ttstub_dsv_t, ttstub_t>      ttstub_ref_t;

    typedef ttstub_ref_t                          reference_t;
    typedef unsigned                              product_t;
    typedef ttstub_dsv_t                          collection_t;

    void setup(const edm::Handle<collection_t>& handle);

    unsigned size() const { return mapping.size(); }

    int get(const reference_t& ref) const;

  private:
    std::map<reference_t, product_t> mapping;
};

#endif
