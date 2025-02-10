#compiler
NVCC       = nvcc
CFLAGS     = -g -G --cudart=shared
TARGET     = fri

#customizable
CUDART_LIB ?= -L/u1/tjandhya/miniforge3/envs/sage/lib
CUDART_INC ?= -Iinclude
LIBS       = -lcudart -lcrypto

#source and objects
SRCS       = params.cu field.cu domain.cu hash-host.cu merkle.cu fft.cu poly.cu poly-eval-launch.cu hash.cu fiat-shamir.cu commit-launch-merkle.cu fri.cu
OBJS       = $(SRCS:.cu=.o)

all: $(TARGET)

#build target
$(TARGET): $(OBJS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(CUDART_LIB) $(LIBS) $(CUDART_INC)

%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

rebuild: clean all