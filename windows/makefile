# Windows build (MSYS2 UCRT64 shell). Run from this folder:
#   make          -> builds windows/moscap_sim.exe
#   make run      -> builds, then runs the simulator from the repo root (writes csv_files/)
#   make plots    -> makes the figures in figs/ from csv_files/
#   make clean    -> removes windows/o_files and windows/moscap_sim.exe
# All Windows build output stays inside this folder, so the Linux
# build in the repo root (o_files/, moscap_sim) is never touched.

# Repo root, relative to this folder
ROOT := ..

# Compiler and flags
CXX := g++
SRC_DIRS := main math poisson carrier ac transient
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra $(addprefix -I$(ROOT)/, $(SRC_DIRS))

# Folders
OBJ_DIR := o_files
BIN := moscap_sim.exe

# Find all .cpp files
SOURCES := $(foreach dir, $(SRC_DIRS), $(wildcard $(ROOT)/$(dir)/*.cpp))

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
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Run from the repo root, since the simulator writes to csv_files/
run: $(BIN)
	cd $(ROOT) && ./windows/$(BIN)

# Make the figures
plots:
	cd $(ROOT)/plot && ./plot.sh

# Clean
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJ_DIR) $(BIN)

# Header dependencies written by the compiler (-MMD), so a changed .hpp rebuilds what includes it
-include $(OBJECTS:.o=.d)

.PHONY: all run plots clean
