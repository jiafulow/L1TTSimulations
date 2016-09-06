#include "DataFormats/Common/interface/Wrapper.h"

#include <map>
#include <vector>
#include <list>
#include <deque>
#include <set>
#include <string>

namespace DataFormats_WrappedStdDictionaries {
  struct dictionary {
  //edm::Wrapper<std::vector<std::vector<int> > >          pippo1;
  edm::Wrapper<std::vector<std::vector<unsigned int> > > pippo2;
  edm::Wrapper<std::vector<std::vector<float> > >        pippo3;
  edm::Wrapper<std::vector<std::vector<bool> > >         pippo4;
  };
}

