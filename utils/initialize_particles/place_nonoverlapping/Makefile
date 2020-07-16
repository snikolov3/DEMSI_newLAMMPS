CXX=g++
CXXFLAGS=-Wall -g #-O2

TARGET = place_particles.out

SOURCES = $(shell echo *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) 

.PHONY: clean

clean:
	$(RM) *.o $(TARGET_EXEC)

