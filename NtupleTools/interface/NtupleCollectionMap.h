#ifndef NtupleTools_NtupleCollectionMap_h_
#define NtupleTools_NtupleCollectionMap_h_

////////////////////////////////////////////////////////////////////////////////
// The main implementation in this file is the CollectionMap class. The other
// classes are there to set up all the types and type deduction methods, so
// that the CollectionMap class works for various input collections.
////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <map>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"


namespace ntuple {

  // Definition for remove_handle
  template <class T>
  struct remove_handle
  {
    typedef T type;
  };

  // Partial specializations for remove_handle
  template <class T>
  struct remove_handle<edm::Handle<T> >
  {
    typedef T type;
  };

  template <class T>
  struct remove_handle<const edm::Handle<T> >
  {
    typedef T type;
  };

  // Collection categories
  struct collection_tag {};
  struct vector_tag : public collection_tag {};
  struct detsetvector_tag : public collection_tag {};
  struct detsetvectornew_tag : public collection_tag {};
  struct trackingparticlecollection_tag : public vector_tag {}; // needs special attention

  // Definition for collection_traits
  template <class Collection>
  struct collection_traits
  {
    typedef typename Collection::collection_category  collection_category;
    typedef typename Collection::collection_type      collection_type;
    typedef typename Collection::identifier_type      identifier_type;
    typedef typename Collection::pointer_type         pointer_type;
    typedef typename Collection::get_identifier       get_identifier;
    typedef typename Collection::get_pointer          get_pointer;
  };

}  // namespace ntuple


// _____________________________________________________________________________
// Collection definition for tracking particles

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

namespace detail {

  struct _trackingparticlecollection {
    typedef ntuple::trackingparticlecollection_tag      collection_category;
    typedef TrackingParticleCollection                  collection_type;
    typedef std::pair<EncodedEventId, unsigned int>     identifier_type;
    typedef TrackingParticleCollection::const_iterator  pointer_type;
    typedef TrackingParticleCollection::const_iterator  iterator1_type;
    typedef TrackingParticle::g4t_iterator              iterator2_type;

    struct get_identifier {  // (geoId, channel)
      identifier_type operator()(iterator1_type it1, iterator2_type it2) {
        return identifier_type(it1->eventId(), it2->trackId());
      }
    };

    struct get_pointer {  // tracking particle
      pointer_type operator()(iterator1_type it1, iterator2_type it2) {
        return pointer_type(it1);
      }
    };
  };

}  // namespace detail


// _____________________________________________________________________________
// Collection definition for tracker digis

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

namespace detail {

  struct _trackerdigicollection {
    typedef Phase2TrackerDigi                         digi_t;

    typedef ntuple::detsetvector_tag                  collection_category;
    typedef edm::DetSetVector<digi_t>                 collection_type;
    typedef std::pair<unsigned int, unsigned int>     identifier_type;
    typedef edm::DetSet<digi_t>::const_iterator       pointer_type;
    typedef edm::DetSetVector<digi_t>::const_iterator iterator1_type;
    typedef edm::DetSet<digi_t>::const_iterator       iterator2_type;

    struct get_identifier {  // (geoId, channel)
      identifier_type operator()(iterator1_type it1, iterator2_type it2) {
        const unsigned channel = Phase2TrackerDigi::pixelToChannel(it2->row(), it2->column());
        return identifier_type(it1->detId(), channel);
      }
    };

    struct get_pointer {  // tracker digi
      pointer_type operator()(iterator1_type it1, iterator2_type it2) {
        return pointer_type(it2);
      }
    };
  };

}  // namespace detail


// _____________________________________________________________________________
// Collection definition for TTClusters, TTStubs

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
//#include "DataFormats/L1TrackTrigger/interface/TTStub.h"

namespace detail {

  struct _ttclustercollection {
    typedef Ref_Phase2TrackerDigi_                            digi_ref_t;
    typedef TTCluster<digi_ref_t>                             ttclus_t;

    typedef ntuple::detsetvector_tag                          collection_category;
    typedef edmNew::DetSetVector<ttclus_t>                    collection_type;
    typedef std::pair<unsigned int, std::vector<digi_ref_t> > identifier_type;
    typedef edmNew::DetSet<ttclus_t>::const_iterator          pointer_type;
    typedef edmNew::DetSetVector<ttclus_t>::const_iterator    iterator1_type;
    typedef edmNew::DetSet<ttclus_t>::const_iterator          iterator2_type;

    struct get_identifier {  // (geoId, hits)
      identifier_type operator()(iterator1_type it1, iterator2_type it2) {
        return identifier_type(it2->getDetId().rawId(), it2->getHits());
      }
    };

    struct get_pointer {  // TTCluster
      pointer_type operator()(iterator1_type it1, iterator2_type it2) {
        return pointer_type(it2);
      }
    };
  };

}  // namespace detail


// _____________________________________________________________________________
// Collection Map
//
// This class loops over an input collection, and builds a map that links
// each element (using a user-defined identifier) to its index in the input
// collection. This index in the input collection is presumably the same index
// the element is stored in the ntuple.
//

template <class Collection>
class CollectionMap {
public:
  // Typedefs
  typedef typename ntuple::collection_traits<Collection>::collection_category collection_category;
  typedef typename ntuple::collection_traits<Collection>::collection_type     collection_type;
  typedef typename ntuple::collection_traits<Collection>::identifier_type     identifier_type;
  typedef typename ntuple::collection_traits<Collection>::pointer_type        pointer_type;
  typedef typename ntuple::collection_traits<Collection>::get_identifier      get_identifier;
  typedef typename ntuple::collection_traits<Collection>::get_pointer         get_pointer;

  typedef edm::Handle<collection_type>                                        handle_type;
  typedef std::map<identifier_type, pointer_type>                             map_type;
  typedef typename map_type::size_type                                        map_size_type;
  typedef int                                                                 index_type;

  // Functions
  void setup(const handle_type& handle);

  pointer_type get(const identifier_type& identifier) const;

  index_type get_index(const identifier_type& identifier) const;

  map_size_type size() const { return mapping.size(); }

private:
  // Details
  void setup_detail(const handle_type& handle, const ntuple::vector_tag&);
  void setup_detail(const handle_type& handle, const ntuple::detsetvector_tag&);
  void setup_detail(const handle_type& handle, const ntuple::detsetvectornew_tag&);
  void setup_detail(const handle_type& handle, const ntuple::trackingparticlecollection_tag&);

  // Member data
  map_type mapping;
};

// _____________________________________________________________________________
template <class Collection>
void CollectionMap<Collection>::setup(const handle_type& handle) {
  setup_detail(handle, collection_category());
}

template <class Collection>
void CollectionMap<Collection>::setup_detail(const handle_type& handle, const ntuple::vector_tag&) {
  if (handle.isValid()) {
    typename ntuple::remove_handle<handle_type>::type::const_iterator it  = handle->begin();
    typename ntuple::remove_handle<handle_type>::type::const_iterator end = handle->end();
    for (; it != end; ++it) {
      const identifier_type& id = get_identifier()(it);
      const pointer_type& pointer = get_pointer()(it);
      mapping.insert(std::make_pair(id, pointer));
    }
  }
}

template <class Collection>
void CollectionMap<Collection>::setup_detail(const handle_type& handle, const ntuple::detsetvector_tag&) {
  if (handle.isValid()) {
    typename ntuple::remove_handle<handle_type>::type::const_iterator it1  = handle->begin();
    typename ntuple::remove_handle<handle_type>::type::const_iterator end1 = handle->end();
    for (; it1 != end1; ++it1) {
      typename ntuple::remove_handle<handle_type>::type::value_type::const_iterator it2  = it1->begin();
      typename ntuple::remove_handle<handle_type>::type::value_type::const_iterator end2 = it1->end();
      for (; it2 != end2; ++it2) {
        const identifier_type& id = get_identifier()(it1, it2);
        const pointer_type& pointer = get_pointer()(it1, it2);
        mapping.insert(std::make_pair(id, pointer));
      }
    }
  }
}

template <class Collection>
void CollectionMap<Collection>::setup_detail(const handle_type& handle, const ntuple::detsetvectornew_tag&) {
  if (handle.isValid()) {
    typename ntuple::remove_handle<handle_type>::type::const_iterator it1  = handle->begin();
    typename ntuple::remove_handle<handle_type>::type::const_iterator end1 = handle->end();
    for (; it1 != end1; ++it1) {
      typename ntuple::remove_handle<handle_type>::type::value_type::const_iterator it2  = it1->begin();
      typename ntuple::remove_handle<handle_type>::type::value_type::const_iterator end2 = it1->end();
      for (; it2 != end2; ++it2) {
        const identifier_type& id = get_identifier()(it1, it2);
        const pointer_type& pointer = get_pointer()(it1, it2);
        mapping.insert(std::make_pair(id, pointer));
      }
    }
  }
}

template <class Collection>
void CollectionMap<Collection>::setup_detail(const handle_type& handle, const ntuple::trackingparticlecollection_tag&) {
  if (handle.isValid()) {
    typename ntuple::remove_handle<handle_type>::type::const_iterator it1  = handle->begin();
    typename ntuple::remove_handle<handle_type>::type::const_iterator end1 = handle->end();
    for (; it1 != end1; ++it1) {
      // Special for trackingparticlecollection_tag: additional loop over g4Tracks
      TrackingParticle::g4t_iterator it2  = it1->g4Track_begin();
      TrackingParticle::g4t_iterator end2 = it1->g4Track_end();
      for (; it2 != end2; ++it2) {
        const identifier_type& id = get_identifier()(it1, it2);
        const pointer_type& pointer = get_pointer()(it1, it2);
        mapping.insert(std::make_pair(id, pointer));
      }
    }
  }
}

template <class Collection>
typename CollectionMap<Collection>::pointer_type CollectionMap<Collection>::get(const identifier_type& identifier) const {
  return mapping.at(identifier);  // must exist
}

template <class Collection>
typename CollectionMap<Collection>::index_type CollectionMap<Collection>::get_index(const identifier_type& identifier) const {
  typename map_type::const_iterator found = mapping.find(identifier);
  if (found == mapping.end())
    return -1;
  index_type index = std::distance(mapping.begin(), found);
  return index;
}

// _____________________________________________________________________________
// Finally

typedef CollectionMap<detail::_trackingparticlecollection>  TrackingParticleCollectionMap;
typedef CollectionMap<detail::_trackerdigicollection>       TrackerDigiCollectionMap;
typedef CollectionMap<detail::_ttclustercollection>         TTClusterCollectionMap;

#endif
