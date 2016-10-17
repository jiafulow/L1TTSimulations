#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "L1TTSimulations/TrackerTools/interface/ModuleIdHelper.h"

using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
using Phase2TrackerTopology    = PixelTopology;


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>

class AnalyzerModuleVertices : public edm::EDAnalyzer {
  public:
    /// Constructor/destructor
    explicit AnalyzerModuleVertices(const edm::ParameterSet&);
    virtual ~AnalyzerModuleVertices();

  private:
    virtual void beginRun(const edm::Run&, const edm::EventSetup&);
    virtual void endRun(const edm::Run&, const edm::EventSetup&);

    virtual void beginJob();
    virtual void endJob();

    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  private:
    /// For event setup
    const TrackerGeometry * theGeometry;
    const TrackerTopology * theTopology;

    /// Maps
    std::map<uint32_t, uint32_t> moduleId0ToStackId;
    std::map<uint32_t, uint32_t> moduleId1ToStackId;
    std::map<uint32_t, uint32_t> moduleId0ToGeoId;
    std::map<uint32_t, uint32_t> moduleId1ToGeoId;

    /// Configurations
    std::string csvfile_;
    int verbose_;
};

AnalyzerModuleVertices::AnalyzerModuleVertices(const edm::ParameterSet& iConfig)
: csvfile_(iConfig.getParameter<std::string>("csv") ),
  verbose_(iConfig.getParameter<int>("verbosity") ) {}

AnalyzerModuleVertices::~AnalyzerModuleVertices() {}

// Here we make the layout of the detector
void AnalyzerModuleVertices::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
    /// Geometry setup
    edm::ESHandle<TrackerGeometry> geometryHandle;
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);
    theGeometry = geometryHandle.product();

    edm::ESHandle<TrackerTopology> topologyHandle;
    iSetup.get<TrackerTopologyRcd>().get(topologyHandle);
    theTopology = topologyHandle.product();

    /// Clear maps
    moduleId0ToStackId.clear();
    moduleId1ToStackId.clear();
    moduleId0ToGeoId.clear();
    moduleId1ToGeoId.clear();

    /// Loop over the detector elements
    for (auto const & det_u : theGeometry->detUnits()) {
        const DetId detId = det_u->geographicalId();

        if (detId.subdetId()!=StripSubdetector::TOB && detId.subdetId()!=StripSubdetector::TID)  // only run on outer tracker
            continue;
        if (!theTopology->isLower(detId))  // loop on the stacks: choose the lower arbitrarily
            continue;

        const DetId stackDetId = theTopology->stack(detId);
        //stackIdToGeoIdMap[stackDetId.rawId()] = detId.rawId();
        //std::cout << theTopology->print(detId) << std::endl;

        const Phase2TrackerGeomDetUnit* pixdet = dynamic_cast<const Phase2TrackerGeomDetUnit*>(det_u);
        assert(pixdet);

        const DetId geoId0 = detId;
        const DetId geoId1 = theTopology->partnerDetId(geoId0);

        uint32_t moduleId0 = ModuleIdHelper::getModuleId(theTopology, detId);
        uint32_t moduleId1 = ModuleIdHelper::getModuleId(theTopology, detId);
        assert(moduleId0 == moduleId1);

        //std::cout << moduleId0 << " " << theTopology->layer(geoId0) << " " << theTopology->tobSide(geoId0) << " " << theTopology->tobRod(geoId0) << " " << theTopology->tobModule(geoId0) << std::endl;
        //std::cout << moduleId0 << " " << theTopology->layer(geoId0) << " " << theTopology->tidSide(geoId0) << " " << theTopology->tidRing(geoId0) << " " << theTopology->tidModule(geoId0) << std::endl;

        if (moduleId0ToGeoId.find(moduleId0) == moduleId0ToGeoId.end() )
            moduleId0ToGeoId.insert(std::make_pair(moduleId0, geoId0.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId0: " << moduleId0 << " geoId0: " << geoId0.rawId() << " existing value in map: " << moduleId0ToGeoId.at(moduleId0) << std::endl;

        if (moduleId1ToGeoId.find(moduleId1) == moduleId1ToGeoId.end() )
            moduleId1ToGeoId.insert(std::make_pair(moduleId1, geoId1.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId1: " << moduleId1 << " geoId1: " << geoId1.rawId() << " existing value in map: " << moduleId1ToGeoId.at(moduleId1) << std::endl;

        if (moduleId0ToStackId.find(moduleId0) == moduleId0ToStackId.end() )
            moduleId0ToStackId.insert(std::make_pair(moduleId0, stackDetId.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId0: " << moduleId0 << " stackId: " << stackDetId.rawId() << " existing value in map: " << moduleId0ToStackId.at(moduleId0) << std::endl;

        if (moduleId1ToStackId.find(moduleId1) == moduleId1ToStackId.end() )
            moduleId1ToStackId.insert(std::make_pair(moduleId1, stackDetId.rawId()) );
        else
            std::cout << "Error: This pair already exists in map! moduleId1: " << moduleId1 << " stackId: " << stackDetId.rawId() << " existing value in map: " << moduleId1ToStackId.at(moduleId1) << std::endl;
    }  // end loop over detector elements
}

void AnalyzerModuleVertices::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {}


void AnalyzerModuleVertices::beginJob() {}


// Here we write into the text file
void AnalyzerModuleVertices::endJob() {
    std::cout << "Map size: " << moduleId0ToStackId.size() << ", " << moduleId1ToStackId.size()
              << ", " << moduleId0ToGeoId.size() << ", " << moduleId1ToGeoId.size() << std::endl;

    // Open text file
    std::ofstream csvfile(csvfile_);

    // Write text file
    csvfile << "moduleId/I, x0_cm/D, y0_cm/D, z0_cm/D, x1_cm/D, y1_cm/D, z1_cm/D, x2_cm/D, y2_cm/D, z2_cm/D, x3_cm/D, y3_cm/D, z3_cm/D" << std::endl;

    for (const auto& kv : moduleId0ToGeoId) {
        const uint32_t moduleId = kv.first;

        for (unsigned i=0; i<2; ++i) {
            const DetId geoId = (i == 0) ? moduleId0ToGeoId.at(kv.first) : moduleId1ToGeoId.at(kv.first);

            const GeomDetUnit* geoUnit = theGeometry->idToDetUnit(geoId);
            const Phase2TrackerGeomDetUnit* pixUnit = dynamic_cast<const Phase2TrackerGeomDetUnit*>(geoUnit);
            const Phase2TrackerTopology* pixTopo = dynamic_cast<const Phase2TrackerTopology*>(&(pixUnit->specificTopology()) );

            const int nrows = pixTopo->nrows();
            const int ncols = pixTopo->ncolumns();

            /// Add 0.5 to get the center of the pixel
            MeasurementPoint mp0(0    , 0    );
            MeasurementPoint mp1(nrows, 0    );
            MeasurementPoint mp2(nrows, ncols);
            MeasurementPoint mp3(0    , ncols);

            /// Find global positions
            const GlobalPoint& pos0 = pixUnit->surface().toGlobal(pixTopo->localPosition(mp0));
            const GlobalPoint& pos1 = pixUnit->surface().toGlobal(pixTopo->localPosition(mp1));
            const GlobalPoint& pos2 = pixUnit->surface().toGlobal(pixTopo->localPosition(mp2));
            const GlobalPoint& pos3 = pixUnit->surface().toGlobal(pixTopo->localPosition(mp3));

            // Positions are in unit of centimeter
            csvfile.unsetf(std::ios_base::floatfield);
            csvfile << moduleId << ", "
                    << std::fixed << std::setprecision(6)
                    << pos0.x() << ", " << pos0.y() << ", " << pos0.z() << ", "
                    << pos1.x() << ", " << pos1.y() << ", " << pos1.z() << ", "
                    << pos2.x() << ", " << pos2.y() << ", " << pos2.z() << ", "
                    << pos3.x() << ", " << pos3.y() << ", " << pos3.z()
                    << std::endl;
        }
    }

    // Close text file
    csvfile.close();

    std::cout << ">>> " << csvfile_ << " is written." << std::endl;
}


// ANALYZE
void AnalyzerModuleVertices::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // Do nothing
}

// DEFINE THIS AS A PLUG-IN
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AnalyzerModuleVertices);
