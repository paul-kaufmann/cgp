appname := cgp

UNAME := $(shell uname)
CXX := /usr/local/bin/g++-11
#CXX := g++-11
CXXFLAGS := -g -Ofast --unroll-loops
LDFLAGS :=
LDLIBS :=
ifeq ($(UNAME), Linux)
	LDLIBS := -lbsd -lm
endif

%srcfiles := $(shell find . -name "*.cpp")
srcfiles := cgp.cpp  cgp-functions.cpp node-function.cpp mutation.cpp optimization_algorithms.cpp statistical_methods.cpp optimizationlog.cpp
objects  := $(patsubst %.cpp, %.o, $(srcfiles))



all: $(appname)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects) $(appname)
	
dist-clean: clean
	rm -f *~ .depend
	
include .depend
