CXX=g++
CXXFLAGS=-c -O3
LDFLAGS=
SOURCES=$(shell find . -name "*.cpp") # searches all subdirectories
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=exe

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o exe -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi90 -lm # link with correct gurobi libraries (...)


.cpp.o:
	$(CXX) $(CXXFLAGS) $<  -o $@ -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi90 -lm

clean:
	rm $(OBJECTS) $(EXECUTABLE)
 
