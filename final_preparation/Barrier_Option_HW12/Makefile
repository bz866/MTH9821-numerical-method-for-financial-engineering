CXXFLAGS += -std=c++1y -Wall -I/usr/include/eigen3 -I/usr/local/include/eigen3 -Wno-error
SOURCES := $(shell find . -name '*.cpp')
HEADERS := $(shell find . -name '*.hpp')
OUTPUT := hw12

all: $(OUTPUT)
$(OUTPUT): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(OUTPUT) $(SOURCES)

clean:
	$(RM) $(OUTPUT) *.o
