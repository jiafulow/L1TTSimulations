#ifndef NtupleTools_NtupleMaker_h_
#define NtupleTools_NtupleMaker_h_

/** \class NtupleMaker
 *
 *  This is stolen from:
 *    https://github.com/BristolTopGroup/NtupleProduction/blob/master/src/RootTupleMakerV2_Tree.cc
 *
 *  Makes a tree out of C++ standard types and vectors of C++ standard types
 *
 *  This class, which is an EDAnalyzer, takes the same "keep" and
 *  "drop" outputCommands parameter as the PoolOutputSource, making a
 *  tree of the selected variables, which it obtains from the EDM
 *  tree.
 *
 *  \author Burt Betchart - University of Rochester <burton.andrew.betchart@cern.ch>
 */

#include "L1TTSimulations/NtupleTools/interface/NtupleCommon.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <map>
#include <string>
#include <vector>
#include "TTree.h"


class NtupleMaker : public edm::EDAnalyzer {
  public:
    // NOTE: STRING is not implemented
    enum LeafType {
        CHAR_T=1, UCHAR_T  , SHORT_T   , USHORT_T  , INT_T  , UINT_T  , INT64_T  , UINT64_T  , LONG64_T  , ULONG64_T  ,
        BOOL_T  , FLOAT_T  , DOUBLE_T  , STRING_T  ,
        CHAR_V  , UCHAR_V  , SHORT_V   , USHORT_V  , INT_V  , UINT_V  , INT64_V  , UINT64_V  , LONG64_V  , ULONG64_V  ,
        BOOL_V  , FLOAT_V  , DOUBLE_V  , STRING_V  ,
        CHAR_V_V, UCHAR_V_V, SHORT_V_V , USHORT_V_V, INT_V_V, UINT_V_V, INT64_V_V, UINT64_V_V, LONG64_V_V, ULONG64_V_V,
        BOOL_V_V, FLOAT_V_V, DOUBLE_V_V, STRING_V_V,
        STRING_INT_M, STRING_UINT_M, STRING_BOOL_M, STRING_FLOAT_M, STRING_DOUBLE_M, STRING_STRING_M,
        NumLeafTypes
    };

    explicit NtupleMaker(const edm::ParameterSet&);

  private:
    //virtual void beginJob() {}
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    //virtual void endJob() {}

    // Connect branch to handle
    class BranchConnector {
      public:
        BranchConnector() {};
        virtual ~BranchConnector() {};
        virtual void connect(const edm::Event&) = 0;
    };

    template <class T>
    class TypedBranchConnector : public BranchConnector {
      public:
        TypedBranchConnector(edm::BranchDescription const*, std::string, TTree*);
        void connect(const edm::Event&);
        void setToken(const edm::EDGetTokenT<T>& token) { token_ = token; }

      private:
        std::string moduleLabel_;
        std::string instanceName_;
        std::string processName_;
        T object_;
        T* ptr_object_;
        edm::EDGetTokenT<T> token_;
    };

    // Register branch
    template <class T>
    void registerBranch(edm::BranchDescription const*, const std::string&);

    // Register branches
    void registerBranches();

    // Member data
    edm::ParameterSet pset;
    TTree * tree;
    std::string treeName;

    std::map<std::string, LeafType> leafmap;
    std::vector<BranchConnector*> connectors;
};

#endif
