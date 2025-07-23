# Makefile for D0 Reconstruction Efficiency Calculator
# 
# This Makefile compiles the D0RecoEfficiency C++ class with std::thread parallelization
# and ROOT libraries for LHCb analysis

# Compiler settings
CXX = /home/niviths/alice/sw/ubuntu2004_x86-64/GCC-Toolchain/v12.2.0-alice1-1/bin/g++
# Try to detect OpenMP support, fallback to no OpenMP if not available
OPENMP_FLAGS := $(shell $(CXX) -fopenmp -E -x c /dev/null >/dev/null 2>&1 && echo "-fopenmp" || echo "")
CXXFLAGS = -Wall -O3 -std=c++17 $(OPENMP_FLAGS)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# Directories
SRCDIR = .
OBJDIR = obj
BINDIR = bin

# Source files
SOURCES = $(SRCDIR)/D0RecoEfficiencyRun.cpp
OBJECTS = $(OBJDIR)/D0RecoEfficiencyRun.o
TARGET = $(BINDIR)/D0RecoEfficiency
TEST_TARGET = $(BINDIR)/test_threading

# Default target
all: directories $(TARGET) $(TEST_TARGET)

# Create necessary directories
directories:
	@mkdir -p $(OBJDIR) $(BINDIR)

# Compile object files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(SRCDIR)/D0RecoEfficiency.h
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -c $< -o $@

# Link executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(ROOTLIBS) -o $@

# Test threading executable
$(TEST_TARGET): $(SRCDIR)/test_threading.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Clean build files
clean:
	rm -rf $(OBJDIR) $(BINDIR)

# Run the program with default parameters
run: $(TARGET)
	./$(TARGET)

# Run with custom parameters (example)
run-custom: $(TARGET)
	./$(TARGET) -i input.root -o output.root -pt 2.0 -eta1 2.0 -eta2 4.5

# Test threading implementation
test: $(TEST_TARGET)
	./$(TEST_TARGET)

# Install dependencies (Ubuntu/Debian)
install-deps:
	sudo apt-get update
	sudo apt-get install build-essential

# Help
help:
	@echo "Available targets:"
	@echo "  all        - Build the program (default)"
	@echo "  clean      - Clean build files"
	@echo "  run        - Run with default parameters"
	@echo "  run-custom - Run with custom parameters"
	@echo "  install-deps - Install dependencies on Ubuntu/Debian"
	@echo "  help       - Show this help"
	@echo "  test       - Test threading implementation"

# Mark phony targets
.PHONY: all clean run run-custom install-deps help directories test

# Compilation flags explanation:
# -Wall: Enable all warnings
# -O3: Optimize for speed
# -std=c++17: Use C++17 standard
# -pthread: Enable POSIX threads (std::thread)
# $(shell root-config --cflags): ROOT compilation flags
# $(shell root-config --libs): ROOT libraries for linking
