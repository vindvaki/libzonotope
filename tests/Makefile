#
# A minimalistic Makefile to get things up and running
# 

#
# Dependencies
# 
LIBS = -lCGAL -lboost_thread -liml -lgmp -lgmpxx
INCLUDE = ../include

#
# Compiler options
#
CXX ?= g++ 
CXXFLAGS = -Wall -Wextra -std=c++11 -O3

#
# Compilation targets
# 
TARGETS = main_cell_enum main_vertex_enum main_zonotope_halfspaces main_zonotope_halfspaces_iml main_zonotope_volume

#
# Rules
#

all: $(TARGETS)

$(TARGETS): %:%.cpp ../include/*.hpp 
	$(CXX) -o $@ $(CXXFLAGS) -I$(INCLUDE) $< $(LIBS)

clean:
	rm -rf $(TARGETS);

.PHONY: all clean
