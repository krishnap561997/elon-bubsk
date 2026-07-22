# ============================================================
# Basilisk Makefile
# Usage:
#   make SRC=filmflow.c TARGET=filmflow
#   make clean
# ============================================================

SHELL := /bin/bash

# ==========================================================
# User Options
# ==========================================================

SRC ?= bubble.c

# BUILD = release (default) or debug
BUILD ?= debug

MODULES = ufrc gcc openmpi imagemagick hdf5

CC99 = gcc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1 -fopenmp

DEBUGFLAGS   = -O2 -g -Wall -disable-dimensions
RELEASEFLAGS = -O3 -fno-strict-aliasing -march=x86-64-v3 -disable-dimensions

# ==========================================================
# Build configuration
# ==========================================================

ifeq ($(BUILD),letitride)
TARGET   = bubble 
QCCFLAGS = $(RELEASEFLAGS)
else
TARGET   = bubble_debug
QCCFLAGS = $(DEBUGFLAGS)
endif

OPENGLIBS = \
	-L$(BASILISK)/gl \
	-lglutils \
	-lfb_tiny

HDFLIBS = \
	-I$$HPC_HDF5_INC \
	-L$$HPC_HDF5_LIB \
	-lhdf5 \
	-lhdf5_hl

LIBS = \
		-I$(BASILISK) \
		$(OPENGLIBS) \
		$(HDFLIBS) \
		-lm

# ==========================================================
# Build
# ==========================================================

all: $(TARGET)

$(TARGET): $(SRC)
	module purge; \
	module load $(MODULES); \
	CC99='$(CC99)' qcc $(QCCFLAGS) \
	$(SRC) -o $(TARGET) \
	$(LIBS)

# ==========================================================
# Run
# ==========================================================

run: $(TARGET)
	sbatch run.sh $(TARGET)

# ==========================================================
# Cleanup
# ==========================================================

clean:
	rm -f filmflow filmflow_debug
	rm -f *.o _*.c

purge:
	rm -f *.dat
	rm -f interface*
	rm -f infc*.dat
	rm -f dump*
	rm -rf images
	rm -rf output
	rm -rf infc 
	rm -rf logs
	rm -rf run_*

.PHONY: all run clean purge
