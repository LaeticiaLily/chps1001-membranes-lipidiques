# Projet CHPS1001 - Analyse de membranes lipidiques
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -fopenmp -Wall -Wextra -std=c++17

all: bin/leaflet lib/libleaflet.so

bin/leaflet: src/leaflet.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

lib/libleaflet.so: src/leaflet.cpp
	@mkdir -p lib
	$(CXX) $(CXXFLAGS) -DLF_LIB -fPIC -shared -o $@ $<

clean:
	rm -rf bin lib
