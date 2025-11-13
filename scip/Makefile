SCIPDIR = $(HOME)/SCIPOptSuite-9.2.4-Linux

CXX = g++
CXXFLAGS = -std=c++17 -O3 -I$(SCIPDIR)/include
LDFLAGS = -L$(SCIPDIR)/lib -lscip -Wl,-rpath,$(SCIPDIR)/lib

TARGET = sc_scip
SOURCES = sc_model_scip.cpp

all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)