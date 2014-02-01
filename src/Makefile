LIBS = -lgmp -lgmpxx
INCLUDE = ../include

CXX ?= g++ 
CXXFLAGS = -Wall -Wextra -std=c++11 -O3

TARGETS = libzonotope_extern.so

all: $(TARGETS)

libzonotope_extern.so: zonotope_volume_extern.cpp ../include/*.hpp
	$(CXX) -shared -rdynamic -fPIC -o $@ $(CXXFLAGS) -I$(INCLUDE) $< $(LIBS)

clean:
	rm -rf $(TARGETS);

.PHONY: all clean
