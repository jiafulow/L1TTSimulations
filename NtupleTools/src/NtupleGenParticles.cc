#include "L1TTSimulations/NtupleTools/interface/NtupleGenParticles.h"


NtupleGenParticles::NtupleGenParticles(const edm::ParameterSet& iConfig) :
  inputTag_(iConfig.getParameter<edm::InputTag>("inputTag")),
  prefix_  (iConfig.getParameter<std::string>("prefix")),
  suffix_  (iConfig.getParameter<std::string>("suffix")),
  selector_(iConfig.existsAs<std::string>("cut") ? iConfig.getParameter<std::string>("cut") : "", true),
  maxN_    (iConfig.getParameter<unsigned>("maxN")) {

    token_ = consumes<reco::GenParticleCollection>(inputTag_);

    produces<std::vector<float> >   (prefix_ + "px"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "py"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "pz"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "E"     + suffix_);
    produces<std::vector<float> >   (prefix_ + "pt"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "eta"   + suffix_);
    produces<std::vector<float> >   (prefix_ + "phi"   + suffix_);
    produces<std::vector<float> >   (prefix_ + "M"     + suffix_);
    produces<std::vector<float> >   (prefix_ + "vx"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "vy"    + suffix_);
    produces<std::vector<float> >   (prefix_ + "vz"    + suffix_);
    produces<std::vector<int16_t> > (prefix_ + "charge"+ suffix_);
    produces<std::vector<int> >     (prefix_ + "pdgId" + suffix_);
    produces<std::vector<int> >     (prefix_ + "status"+ suffix_);
    produces<unsigned>              (prefix_ + "size"  + suffix_);
}

NtupleGenParticles::~NtupleGenParticles() {}

void NtupleGenParticles::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    std::unique_ptr<std::vector<float> >   v_px    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_py    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_pz    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_E     (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_pt    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_eta   (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_phi   (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_M     (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_vx    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_vy    (new std::vector<float>());
    std::unique_ptr<std::vector<float> >   v_vz    (new std::vector<float>());
    std::unique_ptr<std::vector<int16_t> > v_charge(new std::vector<int16_t>());
    std::unique_ptr<std::vector<int> >     v_pdgId (new std::vector<int>());
    std::unique_ptr<std::vector<int> >     v_status(new std::vector<int>());
    std::unique_ptr<unsigned>              v_size  (new unsigned(0));

    //__________________________________________________________________________
    if (!iEvent.isRealData()) {
        edm::Handle<reco::GenParticleCollection> parts;
        //iEvent.getByLabel(inputTag_, parts);
        if (!token_.isUninitialized())
            iEvent.getByToken(token_, parts);

        if (parts.isValid()) {
            edm::LogInfo("NtupleGenParticles") << "Size: " << parts->size();

            unsigned n = 0;
            for (reco::GenParticleCollection::const_iterator it = parts->begin(); it != parts->end(); ++it) {
                if (n >= maxN_)
                    break;
                if (!selector_(*it))
                    continue;

                // Fill the vectors
                v_px    ->push_back(it->px());
                v_py    ->push_back(it->py());
                v_pz    ->push_back(it->pz());
                v_E     ->push_back(it->energy());
                v_pt    ->push_back(it->pt());
                v_eta   ->push_back(it->eta());
                v_phi   ->push_back(it->phi());
                v_M     ->push_back(it->mass());
                v_vx    ->push_back(it->vx());
                v_vy    ->push_back(it->vy());
                v_vz    ->push_back(it->vz());
                v_charge->push_back(it->charge());
                v_pdgId ->push_back(it->pdgId());
                v_status->push_back(it->status());

                n++;
            }
            *v_size = v_px->size();

        } else {
            edm::LogError("NtupleGenParticles") << "Cannot get the product: " << inputTag_;
        }
    }

    //__________________________________________________________________________
    iEvent.put(std::move(v_px)    , prefix_ + "px"    + suffix_);
    iEvent.put(std::move(v_py)    , prefix_ + "py"    + suffix_);
    iEvent.put(std::move(v_pz)    , prefix_ + "pz"    + suffix_);
    iEvent.put(std::move(v_E)     , prefix_ + "E"     + suffix_);
    iEvent.put(std::move(v_pt)    , prefix_ + "pt"    + suffix_);
    iEvent.put(std::move(v_eta)   , prefix_ + "eta"   + suffix_);
    iEvent.put(std::move(v_phi)   , prefix_ + "phi"   + suffix_);
    iEvent.put(std::move(v_M)     , prefix_ + "M"     + suffix_);
    iEvent.put(std::move(v_vx)    , prefix_ + "vx"    + suffix_);
    iEvent.put(std::move(v_vy)    , prefix_ + "vy"    + suffix_);
    iEvent.put(std::move(v_vz)    , prefix_ + "vz"    + suffix_);
    iEvent.put(std::move(v_charge), prefix_ + "charge"+ suffix_);
    iEvent.put(std::move(v_pdgId) , prefix_ + "pdgId" + suffix_);
    iEvent.put(std::move(v_status), prefix_ + "status"+ suffix_);
    iEvent.put(std::move(v_size)  , prefix_ + "size"  + suffix_);
}
