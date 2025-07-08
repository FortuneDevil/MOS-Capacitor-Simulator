# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Icarrier -Imain -Imath -Ipoisson

# Folders
SRC_DIRS := carrier main math poisson
OBJ_DIR := o_files
BIN := moscap_sim

# Find all .cpp files
ALL_SOURCES := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))

# Exclude test files (list here any you don't want to compile)
EXCLUDED_SOURCES := math/testMatrix.cpp math/test.cpp

# Filter out excluded files
SOURCES := $(filter-out $(EXCLUDED_SOURCES), $(ALL_SOURCES))

# Map source files to object files in o_files/
OBJECTS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Default target
all: $(BIN)

# Link executable
$(BIN): $(OBJECTS)
	@echo "Linking: $@"
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into o_files/
$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling: $< -> $@"
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJ_DIR) $(BIN)

.PHONY: all clean
