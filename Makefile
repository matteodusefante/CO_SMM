CURRENT_PATH=.
BENCHMARK=benchmark

# BOOST_LIB=./boost_lib/
# BOOST_BIN=./boost_lib/
# BOOST_P=/usr/local/boost_1_60_0/
# BOOST_B=/usr/local/lib/

SOURCES_BENCHMARK=${BENCHMARK}.cpp
OBJECTS_BENCHMARK=$(SOURCES_BENCHMARK:%.cpp=%.o)

ifeq ($(shell uname),Darwin)
	CXX = g++
else
	CXX = g++
endif

CFLAGS = -Wall -Wextra
CFLAGS += -std=c++11
#CFLAGS += -fopenmp
#CFLAGS += -static -g -S#debug
#CFLAGS += -lpapi
CFLAGS += -O3 -pipe -mfpmath=sse -march=native -funroll-loops -ffast-math
#CFLAGS += -fprefetch-loop-arrays

# LIB = -I $(BOOST_P) #for boost
# BIN = -L $(BOOST_B) #for boost binaries
# UFLAGS = $(BIN) $(LIB)

#TIME = #-lboost_timer
#LIB = -I $(BOOST_P) #for boost
#BIN = -L $(BOOST_B) #for boost binaries
#CFLAGS = -liomp5 -Wall -Wextra -O2 -ffast-math -fopenmp
#CFLAGS += $(BIN) $(LIB)
#CCFLAGS = $(DFLAGS) $(TIME) #$(BIN) #$(LIB)


all: benchmark

benchmark: $(OBJECTS_BENCHMARK)
	$(CXX) -o $(BENCHMARK) $(CFLAGS) $(OBJECTS_BENCHMARK:%=$(CURRENT_PATH)/obj/$(BENCHMARK).o)

$(OBJECTS_BENCHMARK): %.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $(CURRENT_PATH)/obj/$(BENCHMARK).o

.PHONY: clean

clean:
	rm -f obj/*.o
	rm ${BENCHMARK}

