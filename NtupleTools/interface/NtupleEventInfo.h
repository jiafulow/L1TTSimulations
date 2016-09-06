#ifndef NtupleTools_NtupleEventInfo_h_
#define NtupleTools_NtupleEventInfo_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"


class NtupleEventInfo : public edm::EDProducer {
  public:
    explicit NtupleEventInfo(const edm::ParameterSet&);
    ~NtupleEventInfo();

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const std::string   prefix_, suffix_;
};

#endif

