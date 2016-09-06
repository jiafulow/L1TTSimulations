#ifndef NtupleTools_NtupleTTStubs_h_
#define NtupleTools_NtupleTTStubs_h_

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
//#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
//#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
//#include "DataFormats/Common/interface/DetSetVector.h"
//#include "DataFormats/Common/interface/DetSetVectorNew.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"


class NtupleTTStubs : public edm::EDProducer {
  public:
    explicit NtupleTTStubs(const edm::ParameterSet&);
    ~NtupleTTStubs();

    typedef Phase2TrackerDigi                     digi_t;
    typedef PixelDigiSimLink                      digi_simlink_t;
    typedef Ref_Phase2TrackerDigi_                digi_ref_t;

    typedef edm::DetSetVector<digi_t>             digi_dsv_t;
    typedef edm::DetSet<digi_t>                   digi_ds_t;
    typedef edm::DetSetVector<digi_simlink_t>     digi_simlink_dsv_t;
    typedef edm::DetSet<digi_simlink_t>           digi_simlink_ds_t;

    typedef TTCluster<digi_ref_t>                 ttclus_t;
    typedef edmNew::DetSetVector<ttclus_t>        ttclus_dsv_t;
    typedef edmNew::DetSet<ttclus_t>              ttclus_ds_t;
    typedef edm::Ref<ttclus_dsv_t, ttclus_t>      ttclus_ref_t;
    //typedef TTClusterAssociationMap<digi_ref_t>   ttclus_assocmap_t;
    typedef TTStub<digi_ref_t>                    ttstub_t;
    typedef edmNew::DetSetVector<ttstub_t>        ttstub_dsv_t;
    typedef edmNew::DetSet<ttstub_t>              ttstub_ds_t;
    typedef edm::Ref<ttstub_dsv_t, ttstub_t>      ttstub_ref_t;
    //typedef TTStubAssociationMap<digi_ref_t>      ttstub_assocmap_t;

  private:
    //virtual void beginJob();
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob();

    virtual void beginRun(const edm::Run&, const edm::EventSetup&);
    virtual void endRun(const edm::Run&, const edm::EventSetup&);

    // For event setup
    const TrackerGeometry * theGeometry;
    const TrackerTopology * theTopology;
    const MagneticField* theMagneticField;

    const edm::InputTag inputTag_, inputTagMC_, inputTagClus_, inputTagDigi_, inputTagTP_;
    const std::string   prefix_, suffix_;

    StringCutObjectSelector<ttstub_t> selector_;
    const unsigned maxN_;

    edm::EDGetTokenT<ttclus_dsv_t> clusToken_;
    edm::EDGetTokenT<ttstub_dsv_t> stubToken_;
};

#endif
