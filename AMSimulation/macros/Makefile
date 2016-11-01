CXX=`root-config --cxx`
CXXFLAGS=-g -std=c++11 -Wall `root-config --cflags --glibs`

TARGET = analyzeStubs analyzeTracks analyzePatternBank

all: $(TARGET)

% : %.C
	$(CXX) -o $@ $< $(CXXFLAGS)

clean:
	rm -f *.d *.so *.pcm *~
	rm -f $(TARGET)

