# Compiler and flags
NVCC = nvcc
CFLAGS = -g -G -Iinclude
LDFLAGS = -L/u1/tjandhya/miniforge3/envs/sage/lib -lcudart -lcrypto --cudart=shared

# Source and object files
SRC = params.cu field.cu domain.cu hash-host.cu merkle.cu fft.cu poly.cu poly-eval-launch.cu hash.cu fiat-shamir.cu commit-launch-merkle.cu fri.cu
OBJ = $(SRC:.cu=.o)

# Output binary
TARGET = fri

# Paths
SRC_DIR = src
BUILD_DIR = build

# Make sure the build directory exists
$(shell mkdir -p $(BUILD_DIR))

# Pattern rule for compiling .cu files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	$(NVCC) $(CFLAGS) -dc $< -o $@

# Build the executable
$(TARGET): $(addprefix $(BUILD_DIR)/, $(OBJ))
	$(NVCC) -g -G -o $@ $^ $(LDFLAGS)

# Clean rule
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET)

# Default target
all: $(TARGET)

# Run the executable
run: $(TARGET)
	./$(TARGET)