CC=g++
CFLAGS=-c -Wall -O3 \
	   -I/home/mrussell/Applications/boost_1_55_0/include
LDFLAGS=-L/home/mrussell/Applications/boost_1_55_0/lib \
		-lboost_random

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

STATS_SOURCES=stats.cpp
STATS_OBJECTS=$(STATS_SOURCES:.cpp=.o)
STATS_EXECUTABLE=stats

all: $(SOURCES) $(EXECUTABLE) $(STATS_SOURCES) $(STATS_EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

$(STATS_EXECUTABLE): $(STATS_OBJECTS)
	$(CC) $(STATS_OBJECTS) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(STATS_OBJECTS) $(STATS_EXECUTABLE)
