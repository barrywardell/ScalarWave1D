CXXFLAGS=-Wall -g -O2 -I/opt/local/include -fopenmp
LDFLAGS=-L/opt/local/lib -lgsl -lm -lgslcblas
CXX = g++

.PHONY: all

all: wave_bh wave_bh_fo wave_bh_rk2 wave_bh_rk4

wave_bh: wave_bh.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

wave_bh_fo: wave_bh_fo.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

wave_bh_rk2: wave_bh_rk2.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

wave_bh_rk4: wave_bh_rk4.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DFP_PRECISION_DOUBLE -o $@ $<

wave_bh_rk4_extended_double: wave_bh_rk4.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DFP_PRECISION_EXTENDED_DOUBLE -o $@ $<

wave_bh_rk4_quad: wave_bh_rk4.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -DFP_PRECISION_QUAD -lquadmath -o $@ $<

clean:
	rm -rf wave_bh wave_bh_fo wave_bh_rk2 wave_bh_rk4 wave_bh_rk4_extended_double wave_bh_rk4_quad *.dSYM
