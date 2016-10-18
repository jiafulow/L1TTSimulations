#include "L1TTSimulations/NtupleTools/interface/NtupleTrackingVertices.h"


NtupleTrackingVertices::NtupleTrackingVertices(const edm::ParameterSet& iConfig) :
  inputTag_(iConfig.getParameter<edm::InputTag>("inputTag")),
  prefix_  (iConfig.getParameter<std::string>("prefix")),
  suffix_  (iConfig.getParameter<std::string>("suffix")),
  selector_(iConfig.existsAs<std::string>("cut") ? iConfig.getParameter<std::string>("cut") : "", true),
  maxN_    (iConfig.getParameter<unsigned>("maxN")) {

    token_ = consumes<TrackingVertexCollection>(inputTag_);

    produces<std::vector<float> >    (prefix_ + "vx"       + suffix_);
    produces<std::vector<float> >    (prefix_ + "vy"       + suffix_);
    produces<std::vector<float> >    (prefix_ + "vz"       + suffix_);
    produces<std::vector<float> >    (prefix_ + "tof"      + suffix_);
    produces<std::vector<int> >      (prefix_ + "vtxId"    + suffix_);
    produces<std::vector<unsigned> > (prefix_ + "evtId"    + suffix_);
    produces<std::vector<bool> >     (prefix_ + "inVolume" + suffix_);
    produces<std::vector<std::vector<unsigned> > > (prefix_ + "vtxIds" + suffix_);
    produces<unsigned>               (prefix_ + "size"     + suffix_);
}

NtupleTrackingVertices::~NtupleTrackingVertices() {}

void NtupleTrackingVertices::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::unique_ptr<std::vector<float> >    v_vx      (new std::vector<float>());
    std::unique_ptr<std::vector<float> >    v_vy      (new std::vector<float>());
    std::unique_ptr<std::vector<float> >    v_vz      (new std::vector<float>());
    std::unique_ptr<std::vector<float> >    v_tof     (new std::vector<float>());
    std::unique_ptr<std::vector<int> >      v_vtxId   (new std::vector<int>());
    std::unique_ptr<std::vector<unsigned> > v_evtId   (new std::vector<unsigned>());
    std::unique_ptr<std::vector<bool> >     v_inVolume(new std::vector<bool>());
    std::unique_ptr<std::vector<std::vector<unsigned> > > v_vtxIds(new std::vector<std::vector<unsigned> >());
    std::unique_ptr<unsigned>               v_size    (new unsigned(0));

    //__________________________________________________________________________
    if (!iEvent.isRealData()) {
        edm::Handle<TrackingVertexCollection> vertices;
        //iEvent.getByLabel(inputTag_, vertices);
        if (!token_.isUninitialized())
            iEvent.getByToken(token_, vertices);

        if (vertices.isValid()) {
            edm::LogInfo("NtupleTrackingVertices") << "Size: " << vertices->size();

            unsigned n = 0;
            for (TrackingVertexCollection::const_iterator it = vertices->begin(); it != vertices->end(); ++it) {
                if (n >= maxN_)
                    break;
                if (!selector_(*it))
                    continue;

                // Fill the vectors
                const math::XYZTLorentzVectorD& position = it->position();
                v_vx->push_back(position.x());
                v_vy->push_back(position.y());
                v_vz->push_back(position.z());
                v_tof->push_back(position.t());
                int vtxId = it->g4Vertices().empty() ? -99 : it->g4Vertices().begin()->vertexId();
                v_vtxId->push_back(vtxId);
                v_evtId->push_back(it->eventId().rawId());
                v_inVolume->push_back(it->inVolume());

                std::vector<unsigned> vtxIds;  // simVertexId (a.k.a. g4VertexId)
                for (TrackingVertex::g4v_iterator itsim = it->g4Vertices_begin(); itsim != it->g4Vertices_end(); ++itsim) {
                    vtxIds.push_back(itsim->vertexId());
                }
                v_vtxIds->push_back(vtxIds);

                n++;
            }
            *v_size = v_vx->size();

        } else {
            edm::LogError("NtupleTrackingVertices") << "Cannot get the product: " << inputTag_;
        }
    }

    //__________________________________________________________________________
    iEvent.put(std::move(v_vx)      , prefix_ + "vx"       + suffix_);
    iEvent.put(std::move(v_vy)      , prefix_ + "vy"       + suffix_);
    iEvent.put(std::move(v_vz)      , prefix_ + "vz"       + suffix_);
    iEvent.put(std::move(v_tof)     , prefix_ + "tof"      + suffix_);
    iEvent.put(std::move(v_vtxId)   , prefix_ + "vtxId"    + suffix_);
    iEvent.put(std::move(v_evtId)   , prefix_ + "evtId"    + suffix_);
    iEvent.put(std::move(v_inVolume), prefix_ + "inVolume" + suffix_);
    iEvent.put(std::move(v_vtxIds)  , prefix_ + "vtxIds"   + suffix_);
    iEvent.put(std::move(v_size)    , prefix_ + "size"     + suffix_);
}
