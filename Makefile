#PATHS
ifeq ($(PREFIX),)
PREFIX=/usr/local/
endif

CURRENT_DIR 	= $(shell pwd)
LHAPDF      	= $(PREFIX)
BOOST       	= $(PREFIX)
PHOTOSPLINE 	= $(PREFIX)

SOURCES 	= $(wildcard src/*.cpp)

OBJECTS 	= $(SOURCES:.cpp=.o)

INCLUDE_PATH 	= -I/usr/local/include -I./inc
INCLUDE_PATH  	+= -I$(CURRENT_DIR)/inc
INCLUDE_PATH  	+= -I$(CURRENT_DIR)/inc/Dipole_models
INCLUDE_PATH  	+= -I$(LHAPDF)/include
INCLUDE_PATH  	+= -I$(BOOST)/include
INCLUDE_PATH    += -I$(SROOT)/include

#Compiler
# CC 		= clang
# CXX 		= clang++
CC 		= gcc
CXX 		= g++

#Dynamic Library

#Flags
CXX_FLAGS       =  $(INCLUDE_PATH) -I. -O3 -fPIC -std=c++11

# LD 		= clang++
LD 		= g++
LD_FLAGS 	= -L/usr/local/lib/ -L/usr/lib -L$(LHAPDF)/lib -L$(BOOST)/lib -L$(PHOTOSPLINE)/lib
LD_FLAGS    += -L$(SROOT)/lib
LD_FLAGS    += -L$(SROOT)/lib64
LD_FLAGS 	+= -lLHAPDF
LD_FLAGS 	+= -lgsl -lgslcblas
LD_FLAGS	+= -lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex

.PHONY: all clean

CT_OBJ = $(CT)src/CT12Pdf.o $(CURRENT_DIR)src/ct10_xs.o

all: bin/nu_cross.exe bin/nu_cross_classic.exe bin/nu_cross_var.exe bin/nu_total_cross_central.exe bin/nu_cross_full.exe bin/nu_cross_diff.exe
test: bin/test.exe
aaron: bin/nu_cross_full_a_la_aaron_muon.exe bin/nu_cross_full_a_la_aaron_tau.exe

bin/nu_cross.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_classic.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_classic.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_var.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_var.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_total_cross_central.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_total_cross_central.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_full.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_full.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_full_a_la_aaron_muon.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_full_a_la_aaron_muon.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/test.exe: src/lhapdf_cross_section.o src/physconst.o mains/test.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_diff.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_diff.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin/nu_cross_full_a_la_aaron_tau.exe: src/lhapdf_cross_section.o src/physconst.o mains/nu_cross_full_a_la_aaron_tau.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

%.o:%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm src/*.o bin/*.exe
