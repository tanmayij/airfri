CUDA_PATH   ?= /usr/local/cuda
NVCC        := /u1/tjandhya/miniforge3/envs/sage/bin/nvcc
LIBFF_INC   := ../external/libff/
LIBIOP_INC  := ../external/libiop/
LIBFQFFT_INC := ../external/libiop/depends/libfqfft
CXXFLAGS    := -I../additive-fft/C++/ \
	-I../external \
	-O2 -std=c++17 \
	-I../include \
	-I../additive-fft/C++/Cantor \
	-I../additive-fft/C++/libiop/algebra \
	-I$(LIBFF_INC) \
	-I$(LIBFF_INC)/algebra/fields/binary \
	-I$(LIBIOP_INC) \
	-I$(LIBIOP_INC)/algebra/field_subset \
	-I$(LIBFQFFT_INC) \
	-I$(LIBFQFFT_INC)/libfqfft \
	-I$(CONDA_PREFIX)/include
LDFLAGS     := -L../external/libff/build/libff -lff -lssl -lcrypto -lsodium -lgmp -L$(CONDA_PREFIX)/lib
TARGET      := test_fri
SRC         := test_fri.cu ../src/fri_utils.cu ../src/commit_kernel.cu ../src/field.cu ../src/hash.cu ../src/prove.cu ../src/fs_transform.cu ../src/sample_indices.cu \
	../external/libiop/libiop/common/common.cpp \
	../external/libiop/libiop/algebra/utils.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(NVCC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)
