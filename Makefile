#makefile for main, extr, phi, minitest and tinytest


CXX       = clang++
OPTIMIZE  = -fopenmp -std=c++11 -O3 

#OMP_INCL  = -I/usr/local/opt/libomp/include
#LLVM_INCL = -I/usr/local/opt/llvm/include
#FFTW_INCL = -I/usr/local/include

#MRACS used Intel MKL lib for FFT
MKL_ROOT  = /opt/intel/oneapi/mkl/latest
MKL_PATH  = $(MKL_ROOT)/lib
MKL_INCL  = $(MKL_ROOT)/include
MKL_LINK  = -L$(MKL_PATH) -Wl,-rpath,$(MKL_PATH) -lmkl_cdft_core -lmkl_intel_lp64 \
            -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl

#MRACS source files 
SRC_DIR   = src
BUILD_DIR = build
SRCS      = mracs_primary.cpp read_in_simdata.cpp kernel.cpp fourier.cpp
OBJS     := $(subst .cpp,.o,$(SRCS))
SRCS     := $(addprefix $(SRC_DIR)/, $(SRCS))
OBJS     := $(addprefix $(BUILD_DIR)/, $(OBJS))

#MRACS binary files
BINS     := mracs counting

.PHONY: all clean cleanall
all: mracs counting

%: $(OBJS) $(BUILD_DIR)/%.o
	$(CXX) $(OPTIMIZE) $(MKL_LINK) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(OPTIMIZE) -c $^ -o $@ 

cleanall: clean
	rm -f $(BINS)

clean:
	rm -rf $(BUILD_DIR)