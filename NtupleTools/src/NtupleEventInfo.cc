#include "L1TTSimulations/NtupleTools/interface/NtupleEventInfo.h"


NtupleEventInfo::NtupleEventInfo(const edm::ParameterSet& iConfig) :
  prefix_  (iConfig.getParameter<std::string>("prefix")),
  suffix_  (iConfig.getParameter<std::string>("suffix")) {

    produces<unsigned long long> ("event");
    produces<unsigned>           ("run");
    produces<unsigned>           ("lumi");
    produces<int>                ("bx");
    produces<int>                ("orbit");
    produces<unsigned long long> ("time");
}

NtupleEventInfo::~NtupleEventInfo() {}

void NtupleEventInfo::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::unique_ptr<unsigned long long> v_event(new unsigned long long(0));
    std::unique_ptr<unsigned>           v_run  (new unsigned(0));
    std::unique_ptr<unsigned>           v_lumi (new unsigned(0));
    std::unique_ptr<int>                v_bx   (new int(-1));
    std::unique_ptr<int>                v_orbit(new int(-1));
    std::unique_ptr<unsigned long long> v_time (new unsigned long long(0));

    //__________________________________________________________________________
    *v_event = iEvent.id().event();
    *v_run   = iEvent.id().run();
    *v_lumi  = iEvent.id().luminosityBlock();
    *v_bx    = iEvent.bunchCrossing();
    *v_orbit = iEvent.orbitNumber();
    *v_time  = iEvent.time().value();

    //__________________________________________________________________________
    iEvent.put(std::move(v_event), "event");
    iEvent.put(std::move(v_run)  , "run");
    iEvent.put(std::move(v_lumi) , "lumi");
    iEvent.put(std::move(v_bx)   , "bx");
    iEvent.put(std::move(v_orbit), "orbit");
    iEvent.put(std::move(v_time) , "time");
}

