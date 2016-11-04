CXX=`root-config --cxx`

# debug mode
CXXFLAGS=-g -std=c++11 -Wall `root-config --cflags --glibs`

# optimized mode
OPTCXXFLAGS=-O2 -DNDEBUG -std=c++11 -Wall `root-config --cflags --glibs`

TARGET = analyzeStubs analyzeTracks analyzePatternBank runPatternMerging

all: $(TARGET)

%: %.C
	$(CXX) -o $@ $< $(CXXFLAGS)

runPatternMerging: runPatternMerging.C PatternMerging.h PatternMerging.cc
	$(CXX) -o $@ $< $(OPTCXXFLAGS)

clean:
	rm -f *.d *.so *.pcm *~
	rm -f $(TARGET)

