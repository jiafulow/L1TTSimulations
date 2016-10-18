#ifndef NtupleTools_NtupleTrackingVertices_h_
#define NtupleTools_NtupleTrackingVertices_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"


class NtupleTrackingVertices : public edm::EDProducer {
  public:
    explicit NtupleTrackingVertices(const edm::ParameterSet&);
    ~NtupleTrackingVertices();

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const edm::InputTag inputTag_;
    const std::string   prefix_, suffix_;

    StringCutObjectSelector<TrackingVertex> selector_;
    const unsigned maxN_;

    edm::EDGetTokenT<TrackingVertexCollection> token_;
};

#endif
