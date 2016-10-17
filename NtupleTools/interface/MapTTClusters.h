#ifndef NtupleTools_MapTTClusters_h_
#define NtupleTools_MapTTClusters_h_

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include <map>


class MapTTClusters {
  public:
    typedef Ref_Phase2TrackerDigi_                digi_ref_t;
    typedef TTCluster<digi_ref_t>                 ttclus_t;
    typedef edmNew::DetSetVector<ttclus_t>        ttclus_dsv_t;
    typedef edmNew::DetSet<ttclus_t>              ttclus_ds_t;
    typedef edm::Ref<ttclus_dsv_t, ttclus_t>      ttclus_ref_t;

    //typedef ttclus_ref_t                          reference_t;
    typedef std::pair<DetId, std::vector<digi_ref_t> > reference_t;
    typedef unsigned                              product_t;
    typedef ttclus_dsv_t                          collection_t;

    void setup(const edm::Handle<collection_t>& handle);

    unsigned size() const { return mapping.size(); }

    int get(const reference_t& ref) const;

  private:
    std::map<reference_t, product_t> mapping;
};

#endif
