# Windows build (MSYS2 UCRT64 shell). Run from this folder:
#   make          -> builds windows/moscap_sim.exe
#   make run      -> builds, then runs the simulator from the repo root
#   make clean    -> removes windows/o_files and windows/moscap_sim.exe
# All Windows build output stays inside this folder, so the Linux
# build in the repo root (o_files/, moscap_sim) is never touched.

# Repo root, relative to this folder
ROOT := ..

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -I$(ROOT)/carrier -I$(ROOT)/main -I$(ROOT)/math -I$(ROOT)/poisson

# Folders
SRC_DIRS := carrier main math poisson
OBJ_DIR := o_files
BIN := moscap_sim.exe

# Find all .cpp files
ALL_SOURCES := $(foreach dir, $(SRC_DIRS), $(wildcard $(ROOT)/$(dir)/*.cpp))

# Exclude test files (list here any you don't want to compile)
EXCLUDED_SOURCES := $(ROOT)/math/testMatrix.cpp $(ROOT)/math/test.cpp

# Filter out excluded files
SOURCES := $(filter-out $(EXCLUDED_SOURCES), $(ALL_SOURCES))

# Map source files to object files in windows/o_files/
OBJECTS := $(patsubst $(ROOT)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Default target
all: $(BIN)

# Link executable
$(BIN): $(OBJECTS)
	@echo "Linking: $@"
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into windows/o_files/
$(OBJ_DIR)/%.o: $(ROOT)/%.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling: $< -> $@"
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run from the repo root, since the simulator writes to csv_files/
run: $(BIN)
	cd $(ROOT) && ./windows/$(BIN)

# Clean
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJ_DIR) $(BIN)

.PHONY: all run clean
