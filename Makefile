CC=g++
CFLAGS=-c -Wall -Wextra -pedantic -Wno-unused-local-typedefs -O3 -std=c++11 \
	   -I/home/mrussell/Applications/boost_1_55_0/include \
	   -I/tmp/russell-apps/boost_1_55_0/include \
	   -I/usr/include/eigen3
LDFLAGS=-L/home/mrussell/Applications/boost_1_55_0/lib \
		-L/tmp/russell-apps/boost_1_55_0/lib \
		-lboost_program_options

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

ASEP_SOURCES=asep.cpp
ASEP_OBJECTS=$(ASEP_SOURCES:.cpp=.o)
ASEP_EXECUTABLE=asep

STATS_SOURCES=stats.cpp
STATS_OBJECTS=$(STATS_SOURCES:.cpp=.o)
STATS_EXECUTABLE=stats

all: $(SOURCES) $(EXECUTABLE) $(ASEP_SOURCES) $(ASEP_EXECUTABLE) $(STATS_SOURCES) $(STATS_EXECUTABLE)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

$(ASEP_EXECUTABLE): $(ASEP_OBJECTS)
	$(CC) $(ASEP_OBJECTS) $(LDFLAGS) -o $@

$(STATS_EXECUTABLE): $(STATS_OBJECTS)
	$(CC) $(STATS_OBJECTS) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(ASEP_OBJECTS) $(ASEP_EXECUTABLE) $(STATS_OBJECTS) $(STATS_EXECUTABLE)
