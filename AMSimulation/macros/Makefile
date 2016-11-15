CXX=`root-config --cxx`

# debug mode
CXXFLAGS=-g -pthread -pipe -std=c++14 -Wall `root-config --cflags --glibs`

# optimized mode
#OPTCXXFLAGS=-O2 -DNDEBUG -pthread -pipe -std=c++14 -Wall `root-config --cflags --glibs`
OPTCXXFLAGS=-O2 -pthread -pipe -std=c++14 -Wall `root-config --cflags --glibs`

TARGET = analyzeStubs analyzeTracks analyzePatternBank runPatternMerging runRoadMerging

all: $(TARGET)

%: %.C
	$(CXX) -o $@ $< $(CXXFLAGS)

runPatternMerging: runPatternMerging.C PatternMerging.h PatternMerging.cc
	$(CXX) -o $@ $< $(OPTCXXFLAGS)

runRoadMerging: runRoadMerging.C RoadMerging.h RoadMerging.cc
	$(CXX) -o $@ $< $(OPTCXXFLAGS)

clean:
	rm -f *.d *.so *.pcm *~
	rm -f $(TARGET)

