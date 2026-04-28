CXX ?= g++
CC ?= gcc
THIRD_PARTY ?= $(abspath third_party)
WFA2_ROOT ?= $(THIRD_PARTY)/WFA2-lib
ABPOA_ROOT ?= $(THIRD_PARTY)/abPOA
EDLIB_ROOT ?= $(abspath ../longcallD/edlib)

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra
WFA_CPPFLAGS = -I$(WFA2_ROOT)
AB_CPPFLAGS = -I$(ABPOA_ROOT)/include
EDLIB_CPPFLAGS = -I$(EDLIB_ROOT)/include
C_CFLAGS = -O3 -Wall -Wextra
ALIGN_CPPFLAGS = $(WFA_CPPFLAGS) $(AB_CPPFLAGS) $(EDLIB_CPPFLAGS)

WFA2_LIB = $(WFA2_ROOT)/lib/libwfa.a
ABPOA_LIB = $(ABPOA_ROOT)/lib/libabpoa.a

SOURCES_CXX = src/main.cpp \
	src/collect_pipeline.cpp \
	src/bam_digar.cpp \
	src/collect_var.cpp \
	src/collect_phase.cpp \
	src/collect_phase_noisy.cpp \
	src/align.cpp \
	src/collect_output.cpp
SOURCES_C = src/sdust.c src/cgranges.c src/kalloc.c

OBJS = $(SOURCES_CXX:.cpp=.o) $(SOURCES_C:.c=.o)
EDLIB_OBJ = src/edlib.o
OBJS += $(EDLIB_OBJ)

LDFLAGS ?= -lhts -lm -lz -lpthread

.PHONY: all clean check third-party-libs portable-bundle release release-strict

all: pgphase

check: pgphase
	bash scripts/validate_collect_gates.sh

portable-bundle: pgphase
	bash scripts/make_portable_bundle.sh

release: pgphase
	bash scripts/make_release_bundle.sh

release-strict: pgphase
	RUN_CHECKS=1 bash scripts/make_release_bundle.sh

src/align.o: src/align.cpp
	$(CXX) $(CXXFLAGS) $(ALIGN_CPPFLAGS) -c $< -o $@

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

$(EDLIB_OBJ): $(EDLIB_ROOT)/src/edlib.cpp
	$(CXX) $(CXXFLAGS) $(EDLIB_CPPFLAGS) -c $< -o $@

pgphase: $(OBJS) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)

clean:
	rm -f pgphase src/*.o
