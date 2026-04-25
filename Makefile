CXX ?= g++
CC ?= gcc
THIRD_PARTY ?= $(abspath third_party)
WFA2_ROOT ?= $(THIRD_PARTY)/WFA2-lib
ABPOA_ROOT ?= $(THIRD_PARTY)/abPOA

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
WFA_CPPFLAGS = -I$(WFA2_ROOT)
AB_CPPFLAGS = -I$(ABPOA_ROOT)/include
C_CFLAGS = -O3 -Wall -Wextra

WFA2_LIB = $(WFA2_ROOT)/lib/libwfa.a
ABPOA_LIB = $(ABPOA_ROOT)/lib/libabpoa.a

SOURCES_CXX = src/main.cpp src/collect_bam_variation.cpp
SOURCES_C = src/sdust.c src/cgranges.c src/kalloc.c

OBJS = $(SOURCES_CXX:.cpp=.o) $(SOURCES_C:.c=.o)

LDFLAGS ?= -lhts -lm -lz -lpthread

.PHONY: all clean check third-party-libs

all: pgphase

check: pgphase
	bash scripts/validate_collect_gates.sh

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(C_CFLAGS) -c $< -o $@

third-party-libs: $(WFA2_LIB) $(ABPOA_LIB)

$(WFA2_LIB):
	$(MAKE) -C "$(WFA2_ROOT)" setup
	$(MAKE) -C "$(WFA2_ROOT)" lib_wfa

$(ABPOA_LIB):
	$(MAKE) -C "$(ABPOA_ROOT)" libabpoa

pgphase: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f pgphase src/*.o
