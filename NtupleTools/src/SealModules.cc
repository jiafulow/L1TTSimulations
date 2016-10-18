#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "L1TTSimulations/NtupleTools/interface/NtupleEventInfo.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleGenParticles.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleGenEventInfo.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleTrackingParticles.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleTrackingVertices.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleTTStubs.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleMaker.h"

DEFINE_FWK_MODULE(NtupleEventInfo);
DEFINE_FWK_MODULE(NtupleGenParticles);
DEFINE_FWK_MODULE(NtupleGenEventInfo);
DEFINE_FWK_MODULE(NtupleTrackingParticles);
DEFINE_FWK_MODULE(NtupleTrackingVertices);
DEFINE_FWK_MODULE(NtupleTTStubs);
DEFINE_FWK_MODULE(NtupleMaker);

