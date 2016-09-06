#ifndef NtupleTools_NtupleGenParticles_h_
#define NtupleTools_NtupleGenParticles_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


class NtupleGenParticles : public edm::EDProducer {
  public:
    explicit NtupleGenParticles(const edm::ParameterSet&);
    ~NtupleGenParticles();

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const edm::InputTag inputTag_;
    const std::string   prefix_, suffix_;

    StringCutObjectSelector<reco::GenParticle> selector_;
    const unsigned maxN_;

    edm::EDGetTokenT<reco::GenParticleCollection> token_;
};

#endif
