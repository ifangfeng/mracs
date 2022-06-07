#makefile for main, extr, phi, minitest and tinytest

CXX       = clang++
OPTIMIZE  = -fopenmp -std=c++11 -O3 

OMP_PATH  = -L/usr/local/opt/libomp/lib
#LLVM_INCL = -I/usr/local/opt/llvm/include
FFTW_INCL = -I/usr/local/include

#MRACS used Intel MKL lib for FFT
MKL_ROOT  = /opt/intel/oneapi/mkl/latest
MKL_PATH  = $(MKL_ROOT)/lib
MKL_INCL  = $(MKL_ROOT)/include
MKL_LINK  = -L$(MKL_PATH) -Wl,-rpath,$(MKL_PATH) -lmkl_cdft_core -lmkl_intel_lp64 \
			-lmkl_tbb_thread -lmkl_core -ltbb -lstdc++ -lpthread -lm

#MRACS source files 
SRC_DIR   = src
BUILD_DIR = build
SRCS      = mracs_primary.cpp read_in_simdata.cpp kernel.cpp fourier.cpp
OBJS     := $(subst .cpp,.o,$(SRCS))
SRCS     := $(addprefix $(SRC_DIR)/, $(SRCS))
OBJS     := $(addprefix $(BUILD_DIR)/, $(OBJS))

#MRACS binary files
BINS     := mracs counting

.PHONY: all 
.PRECIOUS: $(BUILD_DIR)/%.o
all: mracs counting 

%: $(OBJS) $(BUILD_DIR)/%.o
	$(CXX) $(OPTIMIZE) $(OMP_PATH) $(MKL_LINK) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(OPTIMIZE) $(FFTW_INCL) -c $^ -o $@ 


cleanall: clean
	rm -f $(BINS)

clean:
	rm -f $(BUILD_DIR)/*.o