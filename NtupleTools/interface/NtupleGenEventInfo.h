#ifndef NtupleTools_NtupleGenEventInfo_h_
#define NtupleTools_NtupleGenEventInfo_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


class NtupleGenEventInfo : public edm::EDProducer {
  public:
    explicit NtupleGenEventInfo(const edm::ParameterSet&);
    ~NtupleGenEventInfo();

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    const edm::InputTag genEventInfoTag_;
    const edm::InputTag pileupInfoTag_;
    const edm::InputTag pileupWeightTag_;
    const edm::InputTag pdfWeightTag_;
    const std::string   prefix_, suffix_;

    edm::EDGetTokenT<GenEventInfoProduct>             genEventInfoToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken_;
    edm::EDGetTokenT<double>                          pileupWeightToken_;
    edm::EDGetTokenT<double>                          pdfWeightToken_;
};

#endif

