# Linux build. Run from the repo root:
#   make          -> builds ./moscap_sim
#   make run      -> builds, then runs the simulator (writes csv_files/)
#   make plots    -> makes the figures in figs/ from csv_files/
#   make clean    -> removes o_files/ and moscap_sim
# For Windows (MSYS2), use windows/makefile instead.

# Compiler and flags
CXX := g++
SRC_DIRS := main math poisson carrier ac transient
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra $(addprefix -I, $(SRC_DIRS))

# Folders
OBJ_DIR := o_files
BIN := moscap_sim

# Find all .cpp files
SOURCES := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cpp))

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
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Run the simulator (writes to csv_files/)
run: $(BIN)
	./$(BIN)

# Make the figures
plots:
	cd plot && ./plot.sh

# Clean
clean:
	@echo "Cleaning up..."
	@rm -rf $(OBJ_DIR) $(BIN)

# Header dependencies written by the compiler (-MMD), so a changed .hpp rebuilds what includes it
-include $(OBJECTS:.o=.d)

.PHONY: all run plots clean
