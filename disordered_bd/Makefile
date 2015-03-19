CC=g++
CFLAGS=-c -Wall -Wno-unused-local-typedefs -O3 -std=c++11 \
	   -I/home/mrussell/Applications/boost_1_55_0/include \
	   -I/tmp/russell-apps/boost_1_55_0/include
LDFLAGS=-L/home/mrussell/Applications/boost_1_55_0/lib \
		-L/tmp/russell-apps/boost_1_55_0/lib \
		-lboost_program_options

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
