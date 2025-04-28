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
DICT = LinkDef.cxx
HDRS = LinkDef.h

# Build Rule
all: $(DICT) $(TARGET)

# Dictionary Generation Rule
$(DICT): $(HDRS)
	@echo "Generating ROOT dictionary..."
	rootcling -f $(DICT) -c $(HDRS)

# $(TARGET): $(SRC)
# 	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)


$(TARGET): $(SRC) $(DICT)
	$(CXX) $(CXXFLAGS) $(SRC) $(DICT) -o $@ $(LDFLAGS)

# Clean Rule
clean:
	rm -f $(TARGET) $(DICT) Dict_rdict.pcm
