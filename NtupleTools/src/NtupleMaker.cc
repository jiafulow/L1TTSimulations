#include "L1TTSimulations/NtupleTools/interface/NtupleMaker.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "FWCore/Utilities/interface/TypeID.h"
#include "DataFormats/Provenance/interface/SelectedProducts.h"


template <class T>
NtupleMaker::TypedBranchConnector<T>::TypedBranchConnector(edm::BranchDescription const* desc, std::string type, TTree * tree) :
  moduleLabel_ (desc->moduleLabel() ),
  instanceName_(desc->productInstanceName() ),
  processName_ (desc->processName() ) {

    ptr_object_ = &object_;

    std::string s = instanceName_;
    std::replace(s.begin(), s.end(), '@', '_');
    if (type != "") {
        tree->Branch(s.c_str(), ptr_object_, (s+"/"+type).c_str());  // raw type
    } else {
        tree->Branch(s.c_str(), &ptr_object_);  // STL type
    }
}

template <class T>
void NtupleMaker::TypedBranchConnector<T>::connect(const edm::Event& iEvent) {
    edm::Handle<T> handle;
    //iEvent.getByLabel(moduleLabel_, instanceName_, handle);
    iEvent.getByToken(token_, handle);
    object_ = *handle;
}

template <class T>
void NtupleMaker::registerBranch(edm::BranchDescription const* desc, const std::string& type) {
    TypedBranchConnector<T>* connector = new TypedBranchConnector<T>(desc, type.c_str(), tree);
    connector->setToken(consumes<T>(edm::InputTag(desc->moduleLabel(), desc->productInstanceName())) );
    connectors.push_back(connector);
}

NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) : pset(iConfig) {

    treeName = iConfig.getParameter<std::string>("treeName");

    edm::Service<TFileService> fs;
    tree = fs->make<TTree>(treeName.c_str(), "");

    registerBranches();
}

void NtupleMaker::registerBranches() {

    leafmap[edm::TypeID(typeid(Char_t)).friendlyClassName()] = CHAR_T;
    leafmap[edm::TypeID(typeid(UChar_t)).friendlyClassName()] = UCHAR_T;
    leafmap[edm::TypeID(typeid(Short_t)).friendlyClassName()] = SHORT_T;
    leafmap[edm::TypeID(typeid(UShort_t)).friendlyClassName()] = USHORT_T;
    leafmap[edm::TypeID(typeid(Int_t)).friendlyClassName()] = INT_T;
    leafmap[edm::TypeID(typeid(UInt_t)).friendlyClassName()] = UINT_T;
    leafmap[edm::TypeID(typeid(int64_t)).friendlyClassName()] = INT64_T;
    leafmap[edm::TypeID(typeid(uint64_t)).friendlyClassName()] = UINT64_T;
    leafmap[edm::TypeID(typeid(Long64_t)).friendlyClassName()] = LONG64_T;
    leafmap[edm::TypeID(typeid(ULong64_t)).friendlyClassName()] = ULONG64_T;
    leafmap[edm::TypeID(typeid(Bool_t)).friendlyClassName()] = BOOL_T;
    leafmap[edm::TypeID(typeid(Float_t)).friendlyClassName()] = FLOAT_T;
    leafmap[edm::TypeID(typeid(Double_t)).friendlyClassName()] = DOUBLE_T;
    leafmap[edm::TypeID(typeid(std::string)).friendlyClassName()] = STRING_T;

    leafmap[edm::TypeID(typeid(std::vector<Char_t>)).friendlyClassName()] = CHAR_V;
    leafmap[edm::TypeID(typeid(std::vector<UChar_t>)).friendlyClassName()] = UCHAR_V;
    leafmap[edm::TypeID(typeid(std::vector<Short_t>)).friendlyClassName()] = SHORT_V;
    leafmap[edm::TypeID(typeid(std::vector<UShort_t>)).friendlyClassName()] = USHORT_V;
    leafmap[edm::TypeID(typeid(std::vector<Int_t>)).friendlyClassName()] = INT_V;
    leafmap[edm::TypeID(typeid(std::vector<UInt_t>)).friendlyClassName()] = UINT_V;
    leafmap[edm::TypeID(typeid(std::vector<int64_t>)).friendlyClassName()] = INT64_V;
    leafmap[edm::TypeID(typeid(std::vector<uint64_t>)).friendlyClassName()] = UINT64_V;
    leafmap[edm::TypeID(typeid(std::vector<Long64_t>)).friendlyClassName()] = LONG64_V;
    leafmap[edm::TypeID(typeid(std::vector<ULong64_t>)).friendlyClassName()] = ULONG64_V;
    leafmap[edm::TypeID(typeid(std::vector<Bool_t>)).friendlyClassName()] = BOOL_V;
    leafmap[edm::TypeID(typeid(std::vector<Float_t>)).friendlyClassName()] = FLOAT_V;
    leafmap[edm::TypeID(typeid(std::vector<Double_t>)).friendlyClassName()] = DOUBLE_V;
    leafmap[edm::TypeID(typeid(std::vector<std::string>)).friendlyClassName()] = STRING_V;

    leafmap[edm::TypeID(typeid(std::vector<std::vector<Char_t> >)).friendlyClassName()] = CHAR_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<UChar_t> >)).friendlyClassName()] = UCHAR_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Short_t> >)).friendlyClassName()] = SHORT_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<UShort_t> >)).friendlyClassName()] = USHORT_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Int_t> >)).friendlyClassName()] = INT_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<UInt_t> >)).friendlyClassName()] = UINT_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<int64_t> >)).friendlyClassName()] = INT64_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<uint64_t> >)).friendlyClassName()] = UINT64_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Long64_t> >)).friendlyClassName()] = LONG64_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<ULong64_t> >)).friendlyClassName()] = ULONG64_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Bool_t> >)).friendlyClassName()] = BOOL_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Float_t> >)).friendlyClassName()] = FLOAT_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<Double_t> >)).friendlyClassName()] = DOUBLE_V_V;
    leafmap[edm::TypeID(typeid(std::vector<std::vector<std::string> >)).friendlyClassName()] = STRING_V_V;

    leafmap[edm::TypeID(typeid(std::map<std::string, Int_t>)).friendlyClassName()] = STRING_INT_M;
    leafmap[edm::TypeID(typeid(std::map<std::string, UInt_t>)).friendlyClassName()] = STRING_UINT_M;
    leafmap[edm::TypeID(typeid(std::map<std::string, Bool_t>)).friendlyClassName()] = STRING_BOOL_M;
    leafmap[edm::TypeID(typeid(std::map<std::string, Float_t>)).friendlyClassName()] = STRING_FLOAT_M;
    leafmap[edm::TypeID(typeid(std::map<std::string, Double_t>)).friendlyClassName()] = STRING_DOUBLE_M;
    leafmap[edm::TypeID(typeid(std::map<std::string, std::string>)).friendlyClassName()] = STRING_STRING_M;


    edm::Service<edm::ConstProductRegistry> reg;
    const std::vector<const edm::BranchDescription*>& allBranches = reg->allBranchDescriptions();
    edm::ProductSelectorRules productSelectorRules_(pset, "outputCommands", "NtupleMaker");
    edm::ProductSelector productSelector;
    if (!productSelector.initialized()) {
        productSelector.initialize(productSelectorRules_, allBranches);
    }

    std::set<std::string> branchNames;
    for (std::vector<const edm::BranchDescription*>::const_iterator it = allBranches.begin(); it != allBranches.end(); ++it) {
        const edm::BranchDescription* selection = *it;

        if (productSelector.selected(*selection)) {
            //Check for duplicate branch names
            if (branchNames.find(selection->productInstanceName() ) != branchNames.end() ) {
                throw edm::Exception(edm::errors::Configuration)
                    << "Error in NtupleMaker: More than one branch named: "<< selection->productInstanceName() << std::endl;
            } else {
                branchNames.insert(selection->productInstanceName() );
            }

            // Determine the leaf type
            std::map<std::string, LeafType>::iterator found = leafmap.find(selection->friendlyClassName() );
            switch (found->second) {
                case CHAR_T         : registerBranch<Char_t>                                 (selection, "B"); break;
                case UCHAR_T        : registerBranch<UChar_t>                                (selection, "b"); break;
                case SHORT_T        : registerBranch<Short_t>                                (selection, "S"); break;
                case USHORT_T       : registerBranch<UShort_t>                               (selection, "s"); break;
                case INT_T          : registerBranch<Int_t>                                  (selection, "I"); break;
                case UINT_T         : registerBranch<UInt_t>                                 (selection, "i"); break;
                case INT64_T        : registerBranch<int64_t>                                (selection, "L"); break;
                case UINT64_T       : registerBranch<uint64_t>                               (selection, "l"); break;
                case LONG64_T       : registerBranch<Long64_t>                               (selection, "L"); break;
                case ULONG64_T      : registerBranch<ULong64_t>                              (selection, "l"); break;
                case BOOL_T         : registerBranch<Bool_t>                                 (selection, "O"); break;
                case FLOAT_T        : registerBranch<Float_t>                                (selection, "F"); break;
                case DOUBLE_T       : registerBranch<Double_t>                               (selection, "D"); break;
                case STRING_T       : registerBranch<std::string>                            (selection,  ""); break;
                case CHAR_V         : registerBranch<std::vector<Char_t> >                   (selection,  ""); break;
                case UCHAR_V        : registerBranch<std::vector<UChar_t> >                  (selection,  ""); break;
                case SHORT_V        : registerBranch<std::vector<Short_t> >                  (selection,  ""); break;
                case USHORT_V       : registerBranch<std::vector<UShort_t> >                 (selection,  ""); break;
                case INT_V          : registerBranch<std::vector<Int_t> >                    (selection,  ""); break;
                case UINT_V         : registerBranch<std::vector<UInt_t> >                   (selection,  ""); break;
                case INT64_V        : registerBranch<std::vector<int64_t> >                  (selection,  ""); break;
                case UINT64_V       : registerBranch<std::vector<uint64_t> >                 (selection,  ""); break;
                case LONG64_V       : registerBranch<std::vector<Long64_t> >                 (selection,  ""); break;
                case ULONG64_V      : registerBranch<std::vector<ULong64_t> >                (selection,  ""); break;
                case BOOL_V         : registerBranch<std::vector<Bool_t> >                   (selection,  ""); break;
                case FLOAT_V        : registerBranch<std::vector<Float_t> >                  (selection,  ""); break;
                case DOUBLE_V       : registerBranch<std::vector<Double_t> >                 (selection,  ""); break;
                case STRING_V       : registerBranch<std::vector<std::string> >              (selection,  ""); break;
                case CHAR_V_V       : registerBranch<std::vector<std::vector<Char_t> > >     (selection,  ""); break;
                case UCHAR_V_V      : registerBranch<std::vector<std::vector<UChar_t> > >    (selection,  ""); break;
                case SHORT_V_V      : registerBranch<std::vector<std::vector<Short_t> > >    (selection,  ""); break;
                case USHORT_V_V     : registerBranch<std::vector<std::vector<UShort_t> > >   (selection,  ""); break;
                case INT_V_V        : registerBranch<std::vector<std::vector<Int_t> > >      (selection,  ""); break;
                case UINT_V_V       : registerBranch<std::vector<std::vector<UInt_t> > >     (selection,  ""); break;
                case INT64_V_V      : registerBranch<std::vector<std::vector<int64_t> > >    (selection,  ""); break;
                case UINT64_V_V     : registerBranch<std::vector<std::vector<uint64_t> > >   (selection,  ""); break;
                case LONG64_V_V     : registerBranch<std::vector<std::vector<Long64_t> > >   (selection,  ""); break;
                case ULONG64_V_V    : registerBranch<std::vector<std::vector<ULong64_t> > >  (selection,  ""); break;
                case BOOL_V_V       : registerBranch<std::vector<std::vector<Bool_t> > >     (selection,  ""); break;
                case FLOAT_V_V      : registerBranch<std::vector<std::vector<Float_t> > >    (selection,  ""); break;
                case DOUBLE_V_V     : registerBranch<std::vector<std::vector<Double_t> > >   (selection,  ""); break;
                case STRING_V_V     : registerBranch<std::vector<std::vector<std::string> > >(selection,  ""); break;
                case STRING_INT_M   : registerBranch<std::map<std::string, Int_t> >          (selection,  ""); break;
                case STRING_UINT_M  : registerBranch<std::map<std::string, UInt_t> >         (selection,  ""); break;
                case STRING_BOOL_M  : registerBranch<std::map<std::string, Bool_t> >         (selection,  ""); break;
                case STRING_FLOAT_M : registerBranch<std::map<std::string, Float_t> >        (selection,  ""); break;
                case STRING_DOUBLE_M: registerBranch<std::map<std::string, Double_t> >       (selection,  ""); break;
                case STRING_STRING_M: registerBranch<std::map<std::string, std::string> >    (selection,  ""); break;

                default       : throw edm::Exception(edm::errors::Configuration)
                                    << "Error in NtupleMaker: Cannot handle leaf of type: "
                                    << selection->className()
                                    << " friendlyClassName: "
                                    << selection->friendlyClassName()
                                    << " label: "
                                    << selection->moduleLabel() << ":" << selection->productInstanceName() << ":" << selection->processName()
                                    << std::endl;
            }
        }
    }
}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    for (std::vector<BranchConnector*>::iterator it = connectors.begin(); it != connectors.end(); ++it) {
        (*it)->connect(iEvent);
    }
    tree->Fill();
}
