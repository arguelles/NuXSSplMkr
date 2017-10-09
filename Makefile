#PATHS
CURRENT_DIR 	= $(shell pwd)
TOOLS 	    	= /home/mkroll/Programs/Tools/
LHAPDF      	= $(PREFIX)
BOOST       	= $(PREFIX)
PHOTOSPLINE 	= $(PREFIX)

SOURCES 	= $(wildcard src/*.cpp)

OBJECTS 	= $(SOURCES:.cpp=.o)

INCLUDE_PATH 	= -I/usr/local/include -I./inc
INCLUDE_PATH  	+= -I$(CURRENT_DIR)/inc
INCLUDE_PATH  	+= -I$(CURRENT_DIR)/inc/Dipole_models
INCLUDE_PATH  	+= -I$(TOOLS)/inc
INCLUDE_PATH  	+= -I$(LHAPDF)/include
INCLUDE_PATH  	+= -I$(BOOST)/include
INCLUDE_PATH    += -I$(PHOTOSPLINE)/include

#Libraries
PHOTOSPLINE_LIB     = -L$(PREFIX)/lib
#Compiler
CC 		= clang
CXX 		= clang++

#Dynamic Library

#Flags
CXX_FLAGS       =  $(INCLUDE_PATH) -I. -O3 -fPIC -std=c++11

LD 		= clang++
LD_FLAGS 	= -L/usr/local/lib/ -L/usr/lib -L$(LHAPDF)/lib -L$(BOOST)/lib -L$(TOOLS)/lib
LD_FLAGS	+= $(PHOTOSPLINE_LIB)
LD_FLAGS 	+= -lLHAPDF -lTools
LD_FLAGS 	+= -lgsl -lgslcblas
LD_FLAGS	+= -lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LD_FLAGS	+= -lphotospline

.PHONY: all clean

CT_OBJ = $(CT)src/CT12Pdf.o $(CURRENT_DIR)src/ct10_xs.o

all: bin2/nu_cross.exe bin2/nu_cross_var.exe bin2/nu_cross_simple.exe bin2/nu_cross_simple_hack.exe bin2/nu_cross_full.exe bin2/nu_cross_full_a_la_aaron.exe bin2/nu_cross_full_a_la_aaron_tau.exe
test: bin2/test.exe
diff: bin2/nu_cross_diff.exe

bin2/nu_cross.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_var.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_var.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_simple.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_simple.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_full.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_full.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_full_a_la_aaron.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_full_a_la_aaron.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/test.exe: src/lhapdf_cross_section.o src/physconst.o mains2/test.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_diff.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_diff.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_full_a_la_aaron_tau.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_full_a_la_aaron_tau.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

bin2/nu_cross_simple_hack.exe: src/lhapdf_cross_section.o src/physconst.o mains2/nu_cross_simple_hack.o
	$(LD)  $^ $(LIBS) $(LD_FLAGS) -o $@

%.o:%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm src/*.o
# bin22/*.exe
