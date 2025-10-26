CUDA_PATH   ?= /usr/local/cuda
NVCC        := /u1/tjandhya/miniforge3/envs/sage/bin/nvcc
CXXFLAGS    := -I../include -O2 -std=c++17
TARGET      := libhash.a
SRC         := hash.cu hash_host.cu
OBJ         := hash.o hash_host.o

all: $(TARGET)

%.o: %.cu
	$(NVCC) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	ar rcs $@ $^

clean:
	rm -f $(TARGET) *.o
