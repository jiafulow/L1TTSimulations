#ifndef NtupleTools_NtupleTrackingParticles_h_
#define NtupleTools_NtupleTrackingParticles_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"


class NtupleTrackingParticles : public edm::EDProducer {
  public:
    explicit NtupleTrackingParticles(const edm::ParameterSet&);
    ~NtupleTrackingParticles();

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const edm::InputTag inputTag_;
    const std::string   prefix_, suffix_;

    StringCutObjectSelector<TrackingParticle> selector_;
    const unsigned maxN_;

    edm::EDGetTokenT<TrackingParticleCollection> token_;
};

#endif
