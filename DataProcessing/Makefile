# Compiler
CXX = g++

# ROOT Flags
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)

# Compilation Flags
CXXFLAGS = -Wall -Wextra -O2 $(ROOTCFLAGS)
LDFLAGS  = $(ROOTLIBS)

# Target Executable
TARGET = Thrust.exe
SRC = Thrust.cxx

# Build Rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Clean Rule
clean:
	rm -f $(TARGET)
