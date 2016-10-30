#include "L1TTSimulations/NtupleTools/interface/NtupleTTStubs.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "L1TTSimulations/TrackerTools/interface/ModuleIdHelper.h"
#include "L1TTSimulations/NtupleTools/interface/NtupleCollectionMap.h"

using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
using Phase2TrackerTopology    = PixelTopology;


namespace {

/// findHitGlobalPosition(), findHitGlobalPosition(), findGlobalPosition(), findGlobalDirection() are
/// copied from deprecated Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h
template< typename T >
GlobalPoint findHitGlobalPosition( const Phase2TrackerGeomDetUnit * geomDetUnit, const TTCluster< T > *cluster, unsigned int hitIdx ) {
  /// Add 0.5 to get the center of the pixel
  //const GeomDetUnit* geomDetUnit = idToDetUnit( cluster->getDetId(), cluster->getStackMember() );
  assert(geomDetUnit);
  int row = 0;
  int col = 0;

  if ( cluster->getRows().size() == 0 || cluster->getCols().size() == 0 )
  {
    T hit = cluster->getHits().at(hitIdx);
    row = hit->row();
    col = hit->column();
  }
  else
  {
    row = cluster->getRows().at(hitIdx);
    col = cluster->getCols().at(hitIdx);
  }
  MeasurementPoint mp( row + 0.5, col + 0.5 );
  return geomDetUnit->surface().toGlobal( geomDetUnit->topology().localPosition( mp ) );
}

template<typename T>
GlobalPoint findAverageGlobalPosition( const Phase2TrackerGeomDetUnit * geomDetUnit, const TTCluster< T > *cluster ) {
  double averageX = 0.0;
  double averageY = 0.0;
  double averageZ = 0.0;

  std::vector<T> hits=cluster->getHits();

  /// Loop over the hits and calculate the average coordinates
  if ( hits.size() != 0 )
  {
    for ( unsigned int i = 0; i < hits.size(); i++ )
    {
      GlobalPoint thisHitPosition = findHitGlobalPosition( geomDetUnit, cluster, i );
      averageX += thisHitPosition.x();
      averageY += thisHitPosition.y();
      averageZ += thisHitPosition.z();
    }
    averageX /= hits.size();
    averageY /= hits.size();
    averageZ /= hits.size();
  }
  return GlobalPoint( averageX, averageY, averageZ );
}

template<typename T>
GlobalPoint findGlobalPosition(const Phase2TrackerGeomDetUnit * geomDetUnit, const TTStub<T> *stub) {
  /// Fast version: only inner cluster matters
  return findAverageGlobalPosition( geomDetUnit, stub->getClusterRef(0).get() );
}

template<typename T>
GlobalVector findGlobalDirection(const Phase2TrackerGeomDetUnit * geomDetUnit0, const Phase2TrackerGeomDetUnit * geomDetUnit1, const TTStub<T> *stub) {
  /// Get average position of Clusters composing the Stub
  GlobalPoint innerHitPosition = findAverageGlobalPosition( geomDetUnit0, stub->getClusterRef(0).get() );
  GlobalPoint outerHitPosition = findAverageGlobalPosition( geomDetUnit1, stub->getClusterRef(1).get() );

  /// Calculate the direction
  GlobalVector directionVector( outerHitPosition.x()-innerHitPosition.x(),
                                outerHitPosition.y()-innerHitPosition.y(),
                                outerHitPosition.z()-innerHitPosition.z() );
  return directionVector;
}

static std::map<unsigned, unsigned> stackIdToGeoIdMap;

}  // namespace


NtupleTTStubs::NtupleTTStubs(const edm::ParameterSet& iConfig) :
  inputTag_    (iConfig.getParameter<edm::InputTag>("inputTag")),
  inputTagMC_  (iConfig.getParameter<edm::InputTag>("inputTagMC")),
  inputTagClus_(iConfig.getParameter<edm::InputTag>("inputTagClus")),
  inputTagDigi_(iConfig.getParameter<edm::InputTag>("inputTagDigi")),
  inputTagTP_  (iConfig.getParameter<edm::InputTag>("inputTagTP")),
  prefix_      (iConfig.getParameter<std::string>("prefix")),
  suffix_      (iConfig.getParameter<std::string>("suffix")),
  selector_    (iConfig.existsAs<std::string>("cut") ? iConfig.getParameter<std::string>("cut") : "", true),
  maxN_        (iConfig.getParameter<unsigned>("maxN"))
{
    stubToken_     = consumes<ttstub_dsv_t>              (inputTag_);
    assocmapToken_ = consumes<ttstub_assocmap_t>         (inputTagMC_);
    clusToken_     = consumes<ttclus_dsv_t>              (inputTagClus_);
    digiToken_     = consumes<digi_dsv_t>                (inputTagDigi_);
    digilinkToken_ = consumes<digilink_dsv_t>            (inputTagDigi_);
    tpToken_       = consumes<TrackingParticleCollection>(inputTagTP_);

    // Remember r is rho in cylindrical coordinate system
    produces<std::vector<float> >                   (prefix_ + "x"              + suffix_);
    produces<std::vector<float> >                   (prefix_ + "y"              + suffix_);
    produces<std::vector<float> >                   (prefix_ + "z"              + suffix_);
    produces<std::vector<float> >                   (prefix_ + "r"              + suffix_);
    produces<std::vector<float> >                   (prefix_ + "eta"            + suffix_);
    produces<std::vector<float> >                   (prefix_ + "phi"            + suffix_);
    produces<std::vector<float> >                   (prefix_ + "coordx"         + suffix_);
    produces<std::vector<float> >                   (prefix_ + "coordy"         + suffix_);
    produces<std::vector<float> >                   (prefix_ + "dirx"           + suffix_);
    produces<std::vector<float> >                   (prefix_ + "diry"           + suffix_);
    produces<std::vector<float> >                   (prefix_ + "dirz"           + suffix_);
    produces<std::vector<float> >                   (prefix_ + "roughPt"        + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "modId"          + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "geoId0"         + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "geoId1"         + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "stackId"        + suffix_);
    produces<std::vector<bool> >                    (prefix_ + "barrel"         + suffix_);
    produces<std::vector<bool> >                    (prefix_ + "psmodule"       + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "clusRef0"       + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "clusRef1"       + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "clusWidth0"     + suffix_);
    produces<std::vector<unsigned> >                (prefix_ + "clusWidth1"     + suffix_);
    produces<std::vector<std::vector<unsigned> > >  (prefix_ + "digiChannels"   + suffix_);
    produces<std::vector<std::vector<unsigned> > >  (prefix_ + "digiRefs"       + suffix_);
    produces<std::vector<float> >                   (prefix_ + "trigDist"       + suffix_);
    produces<std::vector<float> >                   (prefix_ + "trigOffset"     + suffix_);
    produces<std::vector<float> >                   (prefix_ + "trigPos"        + suffix_);
    produces<std::vector<float> >                   (prefix_ + "trigBend"       + suffix_);
    produces<std::vector<float> >                   (prefix_ + "separation"     + suffix_);
    produces<std::vector<bool> >                    (prefix_ + "isGenuine"      + suffix_);
    produces<std::vector<bool> >                    (prefix_ + "isUnknown"      + suffix_);
    produces<std::vector<bool> >                    (prefix_ + "isCombinatoric" + suffix_);
    produces<std::vector<int> >                     (prefix_ + "tpId"           + suffix_);
    produces<std::vector<int> >                     (prefix_ + "pdgId"          + suffix_);
    produces<std::vector<float> >                   (prefix_ + "simPt"          + suffix_);
    produces<std::vector<float> >                   (prefix_ + "simEta"         + suffix_);
    produces<std::vector<float> >                   (prefix_ + "simPhi"         + suffix_);
    produces<std::vector<std::vector<int> > >       (prefix_ + "tpIds"          + suffix_);
    //produces<std::vector<std::vector<float> > >     (prefix_ + "fractions"      + suffix_);
    produces<unsigned>                              (prefix_ + "size"           + suffix_);
}

NtupleTTStubs::~NtupleTTStubs() {}

void NtupleTTStubs::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    /// Geometry setup
    edm::ESHandle<TrackerGeometry> geometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
    theGeometry = geometryHandle.product();

    edm::ESHandle<TrackerTopology> topologyHandle;
    iSetup.get<TrackerTopologyRcd>().get(topologyHandle);
    theTopology = topologyHandle.product();

    /// Magnetic field setup
    edm::ESHandle<MagneticField> magneticFieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticFieldHandle);
    theMagneticField = magneticFieldHandle.product();

    // _________________________________________________________________________
    // Make a map for stackId --> geoId
    stackIdToGeoIdMap.clear();
    for (auto const & det_u : theGeometry->detUnits()) {
        DetId detId = det_u->geographicalId();

        if (detId.subdetId()!=StripSubdetector::TOB && detId.subdetId()!=StripSubdetector::TID)  // only run on outer tracker
            continue;
        if (!theTopology->isLower(detId))  // loop on the stacks: choose the lower arbitrarily
            continue;
        DetId stackDetId = theTopology->stack(detId);
        stackIdToGeoIdMap[stackDetId.rawId()] = detId.rawId();
        //std::cout << theTopology->print(detId) << std::endl;

        const Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u);
        assert(pixdet);
    }
}

void NtupleTTStubs::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}

void NtupleTTStubs::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::unique_ptr<std::vector<float> >                  v_x             (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_y             (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_z             (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_r             (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_eta           (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_phi           (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_coordx        (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_coordy        (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_dirx          (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_diry          (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_dirz          (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_roughPt       (new std::vector<float>());
    std::unique_ptr<std::vector<unsigned> >               v_modId         (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_geoId0        (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_geoId1        (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_stackId       (new std::vector<unsigned>());
    std::unique_ptr<std::vector<bool> >                   v_barrel        (new std::vector<bool>());
    std::unique_ptr<std::vector<bool> >                   v_psmodule      (new std::vector<bool>());
    std::unique_ptr<std::vector<unsigned> >               v_clusRef0      (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_clusRef1      (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_clusWidth0    (new std::vector<unsigned>());
    std::unique_ptr<std::vector<unsigned> >               v_clusWidth1    (new std::vector<unsigned>());
    std::unique_ptr<std::vector<std::vector<unsigned> > > v_digiChannels  (new std::vector<std::vector<unsigned> >());
    std::unique_ptr<std::vector<std::vector<unsigned> > > v_digiRefs      (new std::vector<std::vector<unsigned> >());
    std::unique_ptr<std::vector<float> >                  v_trigDist      (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_trigOffset    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_trigPos       (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_trigBend      (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_separation    (new std::vector<float>());
    std::unique_ptr<std::vector<bool> >                   v_isGenuine     (new std::vector<bool>());
    std::unique_ptr<std::vector<bool> >                   v_isUnknown     (new std::vector<bool>());
    std::unique_ptr<std::vector<bool> >                   v_isCombinatoric(new std::vector<bool>());
    std::unique_ptr<std::vector<int> >                    v_tpId          (new std::vector<int>());
    std::unique_ptr<std::vector<int> >                    v_pdgId         (new std::vector<int>());
    std::unique_ptr<std::vector<float> >                  v_simPt         (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_simEta        (new std::vector<float>());
    std::unique_ptr<std::vector<float> >                  v_simPhi        (new std::vector<float>());
    std::unique_ptr<std::vector<std::vector<int> > >      v_tpIds         (new std::vector<std::vector<int> >());
    //std::unique_ptr<std::vector<std::vector<float> > >    v_fractions     (new std::vector<std::vector<float> >());
    std::unique_ptr<unsigned>                             v_size          (new unsigned(0));


    // _________________________________________________________________________
    /// Digi
    edm::Handle<digi_dsv_t> digis;
    //if (inputTagDigi_.encode() != "")
    //    iEvent.getByLabel(inputTagDigi_, digis);
    if (inputTagDigi_.encode() != "" && !digiToken_.isUninitialized()) {
        iEvent.getByToken(digiToken_, digis);
        if (!digis.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTagDigi_;
        }
    }

    TrackerDigiCollectionMap digiMap;
    digiMap.setup(digis);

    /// DigiSimLink
    edm::Handle<digilink_dsv_t> digilinks;
    //if (inputTagDigi_.encode() != "" && !iEvent.isRealData())
    //    iEvent.getByLabel(inputTagDigi_, digilinks);
    if (inputTagDigi_.encode() != "" && !iEvent.isRealData() && !digilinkToken_.isUninitialized()) {
        iEvent.getByToken(digilinkToken_, digilinks);
        if (!digilinks.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTagDigi_;
        }
    }

    /// TrackingParticle
    edm::Handle<TrackingParticleCollection> trackingParticles;
    //if (!iEvent.isRealData())
    //  iEvent.getByLabel(inputTagTP_, trackingParticles);
    if (!iEvent.isRealData() && !tpToken_.isUninitialized()) {
        iEvent.getByToken(tpToken_, trackingParticles);
        if (!trackingParticles.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTagTP_;
        }
    }

    TrackingParticleCollectionMap trkToTPMap;
    trkToTPMap.setup(trackingParticles);

    /// TTCluster
    edm::Handle<ttclus_dsv_t> clusters;
    //iEvent.getByLabel(inputTagClus_, clusters);
    if (!clusToken_.isUninitialized()) {
        iEvent.getByToken(clusToken_, clusters);
        if (!clusters.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTagClus_;
        }
    }

    TTClusterCollectionMap clusMap;
    clusMap.setup(clusters);

    /// TTStub
    edm::Handle<ttstub_dsv_t> stubs;
    //iEvent.getByLabel(inputTag_, stubs);
    if (!stubToken_.isUninitialized()) {
        iEvent.getByToken(stubToken_, stubs);
        if (!digis.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTag_;
        }
    }

    /// MC truth association map
    edm::Handle<ttstub_assocmap_t> assocmap;
    //if (!iEvent.isRealData())
    //    iEvent.getByLabel(inputTagMC_, assocmap);
    if (!iEvent.isRealData() && !assocmapToken_.isUninitialized()) {
        iEvent.getByToken(assocmapToken_, assocmap);
        if (!assocmap.isValid()) {
            edm::LogError("NtupleTTStubs") << "Cannot get the product: " << inputTagMC_;
        }
    }


    if (clusters.isValid() && stubs.isValid()) {
        edm::LogInfo("NtupleTTStubs") << "Size: " << stubs->size();

        unsigned n = 0;
        for (ttstub_dsv_t::const_iterator itv = stubs->begin(); itv != stubs->end(); ++itv) {

            for (ttstub_ds_t::const_iterator it = itv->begin(); it != itv->end(); ++it) {
                //if (n >= maxN_)
                //    break;
                //if (!selector_(*it))
                //    continue;


                const DetId stackDetId = it->getDetId();
                const DetId geoId0     = stackIdToGeoIdMap.at(stackDetId.rawId());
                const DetId geoId1     = theTopology->partnerDetId(geoId0);
                assert(geoId0.det() == DetId::Detector::Tracker);

                uint32_t subdet = geoId0.subdetId();
                bool isBarrel = (subdet == StripSubdetector::TOB);
                bool isEndcap = (subdet == StripSubdetector::TID);
                if (!isBarrel && !isEndcap) {
                    edm::LogError("NtupleTTStubs") << "Inconsistent isBarrel: " << isBarrel << " vs isEndcap: " << isEndcap;
                }
                const TrackerGeometry::ModuleType moduleType = theGeometry->getDetectorType(geoId0);
                bool isPSModule = (moduleType == TrackerGeometry::ModuleType::Ph2PSP) || (moduleType == TrackerGeometry::ModuleType::Ph2PSS);
                bool isSSModule = (moduleType == TrackerGeometry::ModuleType::Ph2SS);
                if (!isPSModule && !isSSModule) {
                    edm::LogError("NtupleTTStubs") << "Inconsistent isPSModule: " << isPSModule << " vs isSSModule: " << isSSModule;
                }
                bool isLower  = theTopology->isLower(geoId0);
                bool isUpper  = theTopology->isUpper(geoId1);
                if (!isLower && !isUpper) {
                    edm::LogError("NtupleTTStubs") << "Inconsistent isLower: " << isLower << " vs isUpper: " << isUpper;
                }

                /// Module ID
                const unsigned moduleId0 = ModuleIdHelper::getModuleId(theTopology, geoId0);
                const unsigned moduleId1 = ModuleIdHelper::getModuleId(theTopology, geoId1);
                edm::LogInfo("NtupleTTStubs") << theTopology->print(geoId0);
                edm::LogInfo("NtupleTTStubs") << "geoId0: " << geoId0.rawId() << " geoId1: " << geoId1.rawId() << " modId0: " << moduleId0 << " modId1: " << moduleId1;
                assert(moduleId0 == moduleId1);

                /// Cluster
                const ttclus_ref_t cluster0 = it->getClusterRef(0);
                const ttclus_ref_t cluster1 = it->getClusterRef(1);
                const TTClusterCollectionMap::identifier_type myClusId0(cluster0->getDetId().rawId(), cluster0->getHits());
                const TTClusterCollectionMap::identifier_type myClusId1(cluster1->getDetId().rawId(), cluster1->getHits());
                const unsigned myClusRef0 = clusMap.get_index(myClusId0);
                const unsigned myClusRef1 = clusMap.get_index(myClusId1);

                /// Topology
                const GeomDetUnit* geoUnit0 = theGeometry->idToDetUnit(geoId0);
                const GeomDetUnit* geoUnit1 = theGeometry->idToDetUnit(geoId1);
                const Phase2TrackerGeomDetUnit* pixUnit0 = dynamic_cast<const Phase2TrackerGeomDetUnit*>(geoUnit0);
                const Phase2TrackerGeomDetUnit* pixUnit1 = dynamic_cast<const Phase2TrackerGeomDetUnit*>(geoUnit1);
                const Surface::PositionType& surfPosition0 = pixUnit0->position();
                const Surface::PositionType& surfPosition1 = pixUnit1->position();
                double separation = (moduleId0 < 110000) ? surfPosition1.perp() - surfPosition0.perp() : ((moduleId0 < 180000) ? surfPosition1.z() - surfPosition0.z() : surfPosition0.z() - surfPosition1.z());

                /// Positions, directions
                const GlobalPoint&      globalPosition  = findGlobalPosition(pixUnit0, &(*it));
                const GlobalVector&     globalDirection = findGlobalDirection(pixUnit0, pixUnit1, &(*it));
                const MeasurementPoint& localCoord0     = cluster0->findAverageLocalCoordinates();
                const MeasurementPoint& localCoord1     = cluster1->findAverageLocalCoordinates();
                //double magfieldStrength = theMagneticField->inTesla(GlobalPoint(0,0,0)).z();
                double roughPt = 0.;  //FIXME: theStackedGeometry->findRoughPt(magfieldStrength, &(*it));
                edm::LogInfo("NtupleTTStubs") << "localCoord0: " << localCoord0.x() << "," << localCoord0.y() << " localCoord1: " << localCoord1.x() << "," << localCoord1.y();

                /// Find cluster widths
                unsigned clusWidth0 = cluster0->findWidth();
                unsigned clusWidth1 = cluster1->findWidth();

                /// Find digis (only the lower cluster)
                const std::vector<int>& theCols = cluster0->getCols();
                const std::vector<int>& theRows = cluster0->getRows();
                assert(theCols.size() != 0 && theRows.size() != 0 && theCols.size() == theRows.size());

                std::vector<unsigned> digiChannels;
                std::vector<unsigned> digiRefs;
                for (unsigned idigi=0; idigi<theCols.size(); ++idigi) {
                    const unsigned channel = Phase2TrackerDigi::pixelToChannel(theRows.at(idigi), theCols.at(idigi));
                    const TrackerDigiCollectionMap::identifier_type map_id(geoId0.rawId(), channel);
                    digiChannels.push_back(channel);
                    digiRefs.push_back(digiMap.get_index(map_id));
                }


                v_x->push_back(globalPosition.x());                   // sviret/HL_LHC: STUB_x
                v_y->push_back(globalPosition.y());                   // sviret/HL_LHC: STUB_y
                v_z->push_back(globalPosition.z());                   // sviret/HL_LHC: STUB_z
                v_r->push_back(globalPosition.perp());
                v_eta->push_back(globalPosition.eta());
                v_phi->push_back(globalPosition.phi());
                v_coordx->push_back(localCoord0.x());                 // sviret/HL_LHC: STUB_strip
                v_coordy->push_back(localCoord0.y());                 // sviret/HL_LHC: STUB_seg
                v_dirx->push_back(globalDirection.x());
                v_diry->push_back(globalDirection.y());
                v_dirz->push_back(globalDirection.z());
                v_roughPt->push_back(roughPt);                        // sviret/HL_LHC: STUB_pt
                v_modId->push_back(moduleId0);
                v_geoId0->push_back(geoId0.rawId());
                v_geoId1->push_back(geoId1.rawId());
                v_stackId->push_back(stackDetId.rawId());
                v_barrel->push_back(isBarrel);
                v_psmodule->push_back(isPSModule);
                v_clusRef0->push_back(myClusRef0);
                v_clusRef1->push_back(myClusRef1);
                v_clusWidth0->push_back(clusWidth0);
                v_clusWidth1->push_back(clusWidth1);
                v_digiChannels->push_back(digiChannels);
                v_digiRefs->push_back(digiRefs);
                v_trigDist->push_back(it->getTriggerDisplacement());  //                            (in full-strip units)
                v_trigOffset->push_back(it->getTriggerOffset());      // sviret/HL_LHC: STUB_cor    (in full-strip units)
                v_trigPos->push_back(it->getTriggerPosition());       //                            (in full-strip units)
                v_trigBend->push_back(it->getTriggerBend());          // sviret/HL_LHC: STUB_deltas (in full-strip units)
                v_separation->push_back(separation);

                // Set to dummy values first
                v_isGenuine->push_back(false);
                v_isUnknown->push_back(false);
                v_isCombinatoric->push_back(false);
                v_tpId->push_back(-1);                                // sviret/HL_LHC: STUB_tp
                v_pdgId->push_back(-99);                              // sviret/HL_LHC: STUB_pdgID
                v_simPt->push_back(-99.);
                v_simEta->push_back(-99.);
                v_simPhi->push_back(-99.);
                v_tpIds->push_back(std::vector<int>());
                //v_fractions->push_back(std::vector<float>());

                /// Retrieve MC association
                if (assocmap.isValid()) {
                    const ttstub_ref_t aref = edmNew::makeRefTo(stubs, it);
                    v_isGenuine->back()      = assocmap->isGenuine(aref);
                    v_isUnknown->back()      = assocmap->isUnknown(aref);
                    v_isCombinatoric->back() = assocmap->isCombinatoric(aref);
                    const edm::Ptr<TrackingParticle> tpptr = assocmap->findTrackingParticlePtr(aref);
                    if (tpptr.isNonnull()) {
                        assert(v_isGenuine->back() == true);
                        v_tpId->back() = tpptr.key();
                        v_pdgId->back() = tpptr->pdgId();
                        v_simPt->back() = tpptr->pt();
                        v_simEta->back() = tpptr->eta();
                        v_simPhi->back() = tpptr->phi();
                    }
                }

                /// Tracking particle association
                if (digis.isValid() && digilinks.isValid()) {
                    std::vector<int> tpIds0;
                    std::vector<int> tpIds1;

                    const std::vector<digi_ref_t>& theHits0 = cluster0->getHits();
                    assert(theHits0.size() == theCols.size());

                    if (digilinks->find(geoId0) != digilinks->end()) {
                        const digilink_ds_t& digilink0 = (*digilinks)[geoId0];
                        for (digilink_ds_t::const_iterator itlink = digilink0.data.begin(); itlink != digilink0.data.end(); ++itlink) {
                            for (unsigned ihit=0; ihit<theHits0.size(); ++ihit) {
                                if (theHits0.at(ihit)->channel() == itlink->channel()) {
                                    if (itlink->fraction() > 0.3) {  // 0.3 is arbitrary
                                        const TrackingParticleCollectionMap::identifier_type map_id(itlink->eventId(), itlink->SimTrackId());
                                        tpIds0.push_back(trkToTPMap.get_index(map_id));
                                    }
                                }
                            }
                        }
                    }

                    const std::vector<digi_ref_t>& theHits1 = cluster1->getHits();

                    if (digilinks->find(geoId1) != digilinks->end()) {
                        const digilink_ds_t& digilink1 = (*digilinks)[geoId1];
                        for (digilink_ds_t::const_iterator itlink = digilink1.data.begin(); itlink != digilink1.data.end(); ++itlink) {
                            for (unsigned ihit=0; ihit<theHits1.size(); ++ihit) {
                                if (theHits1.at(ihit)->channel() == itlink->channel()) {
                                    if (itlink->fraction() > 0.3) {  // 0.3 is arbitrary
                                        const TrackingParticleCollectionMap::identifier_type map_id(itlink->eventId(), itlink->SimTrackId());
                                        tpIds1.push_back(trkToTPMap.get_index(map_id));
                                    }
                                }
                            }
                        }
                    }

                    // Find common elements in two vectors
                    std::vector<int> tpIds;
                    std::sort(tpIds0.begin(), tpIds0.end());
                    std::sort(tpIds1.begin(), tpIds1.end());
                    tpIds0.erase(std::unique(tpIds0.begin(), tpIds0.end()), tpIds0.end());
                    tpIds1.erase(std::unique(tpIds1.begin(), tpIds1.end()), tpIds1.end());
                    std::set_intersection(tpIds0.begin(), tpIds0.end(), tpIds1.begin(), tpIds1.end(), std::back_inserter(tpIds));
                    v_tpIds->back() = tpIds;
                }

                n++;
            }
        }
        *v_size = v_x->size();
    }


    // _________________________________________________________________________
    iEvent.put(std::move(v_x)             , prefix_ + "x"              + suffix_);
    iEvent.put(std::move(v_y)             , prefix_ + "y"              + suffix_);
    iEvent.put(std::move(v_z)             , prefix_ + "z"              + suffix_);
    iEvent.put(std::move(v_r)             , prefix_ + "r"              + suffix_);
    iEvent.put(std::move(v_eta)           , prefix_ + "eta"            + suffix_);
    iEvent.put(std::move(v_phi)           , prefix_ + "phi"            + suffix_);
    iEvent.put(std::move(v_coordx)        , prefix_ + "coordx"         + suffix_);
    iEvent.put(std::move(v_coordy)        , prefix_ + "coordy"         + suffix_);
    iEvent.put(std::move(v_dirx)          , prefix_ + "dirx"           + suffix_);
    iEvent.put(std::move(v_diry)          , prefix_ + "diry"           + suffix_);
    iEvent.put(std::move(v_dirz)          , prefix_ + "dirz"           + suffix_);
    iEvent.put(std::move(v_roughPt)       , prefix_ + "roughPt"        + suffix_);
    iEvent.put(std::move(v_modId)         , prefix_ + "modId"          + suffix_);
    iEvent.put(std::move(v_geoId0)        , prefix_ + "geoId0"         + suffix_);
    iEvent.put(std::move(v_geoId1)        , prefix_ + "geoId1"         + suffix_);
    iEvent.put(std::move(v_stackId)       , prefix_ + "stackId"        + suffix_);
    iEvent.put(std::move(v_barrel)        , prefix_ + "barrel"         + suffix_);
    iEvent.put(std::move(v_psmodule)      , prefix_ + "psmodule"       + suffix_);
    iEvent.put(std::move(v_clusRef0)      , prefix_ + "clusRef0"       + suffix_);
    iEvent.put(std::move(v_clusRef1)      , prefix_ + "clusRef1"       + suffix_);
    iEvent.put(std::move(v_clusWidth0)    , prefix_ + "clusWidth0"     + suffix_);
    iEvent.put(std::move(v_clusWidth1)    , prefix_ + "clusWidth1"     + suffix_);
    iEvent.put(std::move(v_digiChannels)  , prefix_ + "digiChannels"   + suffix_);
    iEvent.put(std::move(v_digiRefs)      , prefix_ + "digiRefs"       + suffix_);
    iEvent.put(std::move(v_trigDist)      , prefix_ + "trigDist"       + suffix_);
    iEvent.put(std::move(v_trigOffset)    , prefix_ + "trigOffset"     + suffix_);
    iEvent.put(std::move(v_trigPos)       , prefix_ + "trigPos"        + suffix_);
    iEvent.put(std::move(v_trigBend)      , prefix_ + "trigBend"       + suffix_);
    iEvent.put(std::move(v_separation)    , prefix_ + "separation"     + suffix_);
    iEvent.put(std::move(v_isGenuine)     , prefix_ + "isGenuine"      + suffix_);
    iEvent.put(std::move(v_isUnknown)     , prefix_ + "isUnknown"      + suffix_);
    iEvent.put(std::move(v_isCombinatoric), prefix_ + "isCombinatoric" + suffix_);
    iEvent.put(std::move(v_tpId)          , prefix_ + "tpId"           + suffix_);
    iEvent.put(std::move(v_pdgId)         , prefix_ + "pdgId"          + suffix_);
    iEvent.put(std::move(v_simPt)         , prefix_ + "simPt"          + suffix_);
    iEvent.put(std::move(v_simEta)        , prefix_ + "simEta"         + suffix_);
    iEvent.put(std::move(v_simPhi)        , prefix_ + "simPhi"         + suffix_);
    iEvent.put(std::move(v_tpIds)         , prefix_ + "tpIds"          + suffix_);
    //iEvent.put(std::move(v_fractions)     , prefix_ + "fractions"      + suffix_);
    iEvent.put(std::move(v_size)          , prefix_ + "size"           + suffix_);
}
