#makefile for main, extr, phi, minitest and tinytest
CXX       = clang++
CXXFLAGES = -std=c++11 -O3
FOMP      = -fopenmp
OMP_INCL  = -I/usr/local/opt/libomp/include
LLVM_INCL = -I/usr/local/opt/llvm/include
FFTW_INCL = -I/usr/local/include
OBJS      = mracs_primary.o read_in_simdata.o kernel.o fourier.o
MKL_LINK  = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm

#SRCS=tool.cc support.cc
#OBJS=$(subst .cc,.o,$(SRCS))
all: mracs main counting
.PHONY: all

mracs: $(OBJS) mracs.o
	$(CXX) $(FOMP) $(CXXFLAGES) $(OBJS) mracs.o -o mracs $(LLVM_INCL) $(OMP_INCL) $(FFTW_INCL) $(MKL_LINK)
main: main.o kernel.o fourier.o mracs_primary.o
	$(CXX) $(CXXFLAGES) main.o kernel.o fourier.o mracs_primary.o -o main
counting: $(OBJS) counting.o
	$(CXX) $(FOMP) $(CXXFLAGES) $(OBJS) counting.o -o counting $(LLVM_INCL) $(OMP_INCL) $(FFTW_INCL) $(MKL_LINK)


mracs.o: mracs.cpp mracs.hpp
	$(CXX) $(FOMP) $(CXXFLAGES) -c mracs.cpp $(LLVM_INCL) $(OMP_INCL) $(FFTW_INCL)
mracs_primary.o: mracs_primary.cpp mracs.hpp
	$(CXX) $(FOMP) $(CXXFLAGES) -c mracs_primary.cpp $(LLVM_INCL) $(OMP_INCL) $(FFTW_INCL)
read_in_simdata.o: read_in_simdata.cpp mracs.hpp
	$(CXX) $(CXXFLAGES) -c read_in_simdata.cpp $(FFTW_INCL)
main.o: main.cpp mracs.hpp
	$(CXX) $(CXXFLAGES) -c main.cpp
counting.o: counting.cpp mracs.hpp
	$(CXX) $(FOMP) $(CXXFLAGES) -c counting.cpp $(LLVM_INCL) $(OMP_INCL) $(FFTW_INCL)
kernel.o: kernel.cpp
	$(CXX) $(CXXFLAGES) -c kernel.cpp
fourier.o: fourier.cpp
	$(CXX) $(CXXFLAGES) -c fourier.cpp


.PHONY: cleanall clean 

cleanall: clean
	rm -f mracs main counting
clean:
	rm -f *.o