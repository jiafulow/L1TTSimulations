#include "L1TTSimulations/NtupleTools/interface/NtupleGenEventInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


NtupleGenEventInfo::NtupleGenEventInfo(const edm::ParameterSet& iConfig) :
  genEventInfoTag_(iConfig.getParameter<edm::InputTag>("genEventInfo")),
  pileupInfoTag_(iConfig.getParameter<edm::InputTag>("pileupInfo")),
  pileupWeightTag_(iConfig.getParameter<edm::InputTag>("pileupWeight")),
  pdfWeightTag_(iConfig.getParameter<edm::InputTag>("pdfWeight")),
  prefix_  (iConfig.getParameter<std::string>("prefix")),
  suffix_  (iConfig.getParameter<std::string>("suffix")) {

    genEventInfoToken_ = consumes<GenEventInfoProduct>(genEventInfoTag_);
    pileupInfoToken_   = consumes<std::vector<PileupSummaryInfo> >(pileupInfoTag_);
    pileupWeightToken_ = mayConsume<double>(pileupWeightTag_);
    pdfWeightToken_    = mayConsume<double>(pdfWeightTag_);

    produces<int>                     (prefix_ + "nPV"           + suffix_);
    produces<std::vector<int> >       (prefix_ + "puNtrks"       + suffix_);
    produces<std::vector<float> >     (prefix_ + "puZpositions"  + suffix_);
    produces<float>                   (prefix_ + "trueNPV"       + suffix_);
    produces<float>                   (prefix_ + "weightMC"      + suffix_);
    produces<float>                   (prefix_ + "weightPU"      + suffix_);
    produces<float>                   (prefix_ + "weightPDF"     + suffix_);
}

NtupleGenEventInfo::~NtupleGenEventInfo() {}

void NtupleGenEventInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::unique_ptr<int>                    v_nPV          (new int(-99.));
    std::unique_ptr<std::vector<int> >      v_puNtrks      (new std::vector<int>());
    std::unique_ptr<std::vector<float> >    v_puZpositions (new std::vector<float>());
    std::unique_ptr<float>                  v_trueNPV      (new float(-99.));
    std::unique_ptr<float>                  v_weightMC     (new float(-99.));
    std::unique_ptr<float>                  v_weightPU     (new float(-99.));
    std::unique_ptr<float>                  v_weightPDF    (new float(-99.));

    //__________________________________________________________________________
    if (!iEvent.isRealData()) {
        edm::Handle<GenEventInfoProduct> genEventInfo;
        //iEvent.getByLabel(genEventInfoTag_, genEventInfo);
        if (!genEventInfoToken_.isUninitialized())
            iEvent.getByToken(genEventInfoToken_, genEventInfo);
        if (genEventInfo.isValid()) {
            *v_weightMC = genEventInfo->weight();
        }

        edm::Handle<std::vector<PileupSummaryInfo> > pileupInfo;
        //iEvent.getByLabel(pileupInfoTag_, pileupInfo);
        if (!pileupInfoToken_.isUninitialized())
            iEvent.getByToken(pileupInfoToken_, pileupInfo);
        if (pileupInfo.isValid()) {
            for(std::vector<PileupSummaryInfo>::const_iterator it = pileupInfo->begin(); it != pileupInfo->end(); ++it) {
                if (it->getBunchCrossing()==0) {
                    *v_nPV          = it->getPU_NumInteractions();
                    *v_puNtrks      = it->getPU_ntrks_highpT();
                    *v_puZpositions = it->getPU_zpositions();
                    *v_trueNPV      = it->getTrueNumInteractions();
                    break;
                }
            }
        }

        if (pileupWeightTag_.encode() != "") {
            edm::Handle<double> pileupWeight;
            //iEvent.getByLabel(pileupWeightTag_, pileupWeight);
            if (!pileupWeightToken_.isUninitialized())
                iEvent.getByToken(pileupWeightToken_, pileupWeight);
            if (pileupWeight.isValid()) {
                *v_weightPU = *pileupWeight;
            }
        }

        if (pdfWeightTag_.encode() != "") {
            edm::Handle<double> pdfWeight;
            //iEvent.getByLabel(pdfWeightTag_, pdfWeight);
            if (!pdfWeightToken_.isUninitialized())
                iEvent.getByToken(pdfWeightToken_, pdfWeight);
            if (pdfWeight.isValid()) {
                *v_weightPDF = *pdfWeight;
            }
        }
    }

    //__________________________________________________________________________
    iEvent.put(std::move(v_nPV)         , prefix_ + "nPV"          + suffix_);
    iEvent.put(std::move(v_puNtrks)     , prefix_ + "puNtrks"      + suffix_);
    iEvent.put(std::move(v_puZpositions), prefix_ + "puZpositions" + suffix_);
    iEvent.put(std::move(v_trueNPV)     , prefix_ + "trueNPV"      + suffix_);
    iEvent.put(std::move(v_weightMC)    , prefix_ + "weightMC"     + suffix_);
    iEvent.put(std::move(v_weightPU)    , prefix_ + "weightPU"     + suffix_);
    iEvent.put(std::move(v_weightPDF)   , prefix_ + "weightPDF"    + suffix_);
}

