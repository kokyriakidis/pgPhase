CXX ?= g++
CC ?= gcc
THIRD_PARTY ?= $(abspath third_party)
WFA2_ROOT ?= $(THIRD_PARTY)/WFA2-lib
ABPOA_ROOT ?= $(THIRD_PARTY)/abPOA
EDLIB_ROOT ?= $(THIRD_PARTY)/edlib/edlib
GBZ_BASE_ROOT ?= $(THIRD_PARTY)/gbz-base
# Prefer the rustup-managed toolchain over any system-installed rustc/cargo.
CARGO       ?= $(firstword $(wildcard $(HOME)/.cargo/bin/cargo) cargo)
CARGO_RUSTC ?= $(firstword $(wildcard $(HOME)/.cargo/bin/rustc) rustc)

GBZ_QUERY_BIN = $(GBZ_BASE_ROOT)/target/release/query
GAF2DB_BIN    = $(GBZ_BASE_ROOT)/target/release/gaf2db
GBZ2DB_BIN    = $(GBZ_BASE_ROOT)/target/release/gbz2db

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -MMD -MP
WFA_CPPFLAGS = -I$(WFA2_ROOT)
AB_CPPFLAGS = -I$(ABPOA_ROOT)/include
EDLIB_CPPFLAGS = -I$(EDLIB_ROOT)/include
C_CFLAGS = -O3 -Wall -Wextra
ALIGN_CPPFLAGS = $(WFA_CPPFLAGS) $(AB_CPPFLAGS) $(EDLIB_CPPFLAGS)

WFA2_LIB = $(WFA2_ROOT)/lib/libwfa.a
ABPOA_LIB = $(ABPOA_ROOT)/lib/libabpoa.a

GBZ_FFI_DIR = src/gbz_ffi
GBZ_FFI_LIB = $(GBZ_FFI_DIR)/target/release/libpgphase_gbz_ffi.a
# System libraries required by the Rust static library (SQLite is bundled).
GBZ_FFI_SYSLIBS = -ldl -lrt

SOURCES_CXX = src/main.cpp \
	src/collect_pipeline.cpp \
	src/region_excludes.cpp \
	src/build_catalog.cpp \
	src/graph_collect.cpp \
	src/hybrid_collect.cpp \
	src/hybrid_inject.cpp \
	src/bam_digar.cpp \
	src/noise_filter.cpp \
	src/collect_var.cpp \
	src/collect_phase.cpp \
	src/collect_phase_pgbam.cpp \
	src/collect_phase_noisy.cpp \
	src/align.cpp \
	src/collect_bam_output.cpp \
	src/collect_output.cpp \
	src/graph_bam_adapter.cpp \
	src/graph_sites.cpp \
	src/graph_query.cpp
SOURCES_C = src/sdust.c src/cgranges.c src/kalloc.c

OBJS = $(SOURCES_CXX:.cpp=.o) $(SOURCES_C:.c=.o)
EDLIB_OBJ = src/edlib.o
OBJS += $(EDLIB_OBJ)

LDFLAGS ?= -lhts -lm -lz -lpthread

-include $(patsubst %.cpp,%.d,$(SOURCES_CXX))

.PHONY: all clean check unit-tests third-party-libs gbz-base portable-bundle release release-strict

all: pgphase

check: pgphase
	bash scripts/validate_collect_gates.sh

unit-tests: test_phase_block_stitch test_graph_sites test_graph_bam_adapter test_hybrid_inject test_noise_filter test_region_exclude
	./test_phase_block_stitch
	./test_graph_sites
	./test_graph_bam_adapter
	./test_hybrid_inject
	./test_noise_filter
	./test_region_exclude

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

gbz-base: $(GBZ_QUERY_BIN) $(GAF2DB_BIN)

$(GBZ_QUERY_BIN) $(GAF2DB_BIN) $(GBZ2DB_BIN):
	cd $(GBZ_BASE_ROOT) && RUSTC=$(CARGO_RUSTC) $(CARGO) build --release --bin query --bin gaf2db --bin gbz2db

$(WFA2_LIB):
	$(MAKE) -C "$(WFA2_ROOT)" setup
	$(MAKE) -C "$(WFA2_ROOT)" lib_wfa

$(ABPOA_LIB):
	$(MAKE) -C "$(ABPOA_ROOT)" libabpoa

$(EDLIB_OBJ): $(EDLIB_ROOT)/src/edlib.cpp
	$(CXX) $(CXXFLAGS) $(EDLIB_CPPFLAGS) -c $< -o $@

$(GBZ_FFI_LIB): $(wildcard $(GBZ_FFI_DIR)/lib.rs $(GBZ_FFI_DIR)/Cargo.toml)
	cd $(GBZ_FFI_DIR) && RUSTC=$(CARGO_RUSTC) $(CARGO) build --release

pgphase: $(OBJS) $(WFA2_LIB) $(ABPOA_LIB) $(GBZ_FFI_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(WFA2_LIB) $(ABPOA_LIB) $(GBZ_FFI_LIB) $(LDFLAGS) $(GBZ_FFI_SYSLIBS)

test_phase_block_stitch: src/test_phase_block_stitch.cpp src/collect_phase.o src/collect_phase_pgbam.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< src/collect_phase.o src/collect_phase_pgbam.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)

test_graph_sites: src/test_graph_sites.cpp src/graph_sites.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_region_exclude: src/test_region_exclude.cpp src/region_excludes.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_graph_bam_adapter: src/test_graph_bam_adapter.cpp src/graph_bam_adapter.o src/noise_filter.o src/graph_sites.o src/graph_query.o src/collect_phase.o src/collect_phase_pgbam.o src/collect_output.o $(GBZ_FFI_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ src/cgranges.o src/kalloc.o src/sdust.o $(LDFLAGS) $(GBZ_FFI_SYSLIBS)

test_noise_filter: src/test_noise_filter.cpp src/graph_bam_adapter.o src/noise_filter.o src/graph_sites.o src/graph_query.o src/collect_phase.o src/collect_phase_pgbam.o src/collect_output.o $(GBZ_FFI_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^ src/cgranges.o src/kalloc.o src/sdust.o $(LDFLAGS) $(GBZ_FFI_SYSLIBS)

test_hybrid_inject: src/test_hybrid_inject.cpp src/hybrid_inject.o src/collect_phase.o src/collect_phase_pgbam.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $< src/hybrid_inject.o src/collect_phase.o src/collect_phase_pgbam.o src/collect_phase_noisy.o src/collect_output.o src/collect_var.o src/noise_filter.o src/align.o src/cgranges.o src/kalloc.o src/sdust.o $(EDLIB_OBJ) $(WFA2_LIB) $(ABPOA_LIB) $(LDFLAGS)

clean:
	rm -f pgphase test_phase_block_stitch test_graph_sites test_graph_bam_adapter test_hybrid_inject test_noise_filter test_region_exclude src/*.o src/*.d
