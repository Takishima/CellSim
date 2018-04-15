# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Main makefile for CellSim
# Nguyen Damien
# October 2012 - May 2013
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# shell environment
export SHELL = /bin/bash

SED = gsed
UNAME := $(shell uname)

# Use this variable to add path to your dynamic libraries
export LDFLAGS=

# =============================================================================
# Choose compiler

CXX = g++
CXX = clang++
#CXX = cc_args.py clang++

# =============================================================================
# Add relevant library paths

LDFLAGS = -L/opt/local/lib

# ==============================================================================
# Compilation defines

DEFINES += WITH_RANGE_CHECK
DEFINES += MULTI_THREAD
DEFINES += HALIDI_2009
#DEFINES += HALIDI_2012
DEFINES += NO_UNITS

# ==============================================================================
# Below this line, everything is automatic !
# =============================================================================

# Programs to create
# BIFDIAG_EXEC = $(BUILDDIR)/bifdiag.x
# HALIDI_EXEC = $(BUILDDIR)/halidi.x
# KOENIGS_EXEC = $(BUILDDIR)/koenigs.x
MAIN_EXEC = $(BUILDDIR)/cellsim.x
export TEST_EXEC = $(BUILDDIR)/test.x

# Subdirectories
BUILDDIR = bin
COREDIR = core
DOCDIR = documentation
FILEDIR = file
INTEGRATORDIR = integrators
LIBDIR = libs
MODELDIR = model
PROGRAMSDIR = programs
TESTDIR = test

make_bin := $(shell test -e $(BUILDDIR) || mkdir -p $(BUILDDIR))

# name of subdirectory for objects in all modules
export OBJDIR := $(shell echo `git rev-parse --abbrev-ref HEAD`-branch-obj)

#=======================================
# Figure out which compiler we are using

CXX_VERSION := $(shell $(CXX) --version)
export IS_GCC := $(shell $(CXX) --version | grep gcc > /dev/null; echo $$?)
export IS_CLANG := $(shell $(CXX) --version | grep clang > /dev/null; echo $$?)

#=======================================
# Get objects list from sub-makefiles
SHELL_EXPORT = "OBJDIR=$(OBJDIR)"

COREOBJ = $(shell make -sC $(COREDIR) -f Make.obj getobj $(SHELL_EXPORT))
INTEGRATOROBJ = $(shell make -sC $(INTEGRATORDIR) -f Make.obj getobj $(SHELL_EXPORT))
FILEOBJ = $(shell make -sC $(FILEDIR) -f Make.obj getobj $(SHELL_EXPORT))
MODELOBJ = $(shell make -sC $(MODELDIR) -f Make.obj getobj $(SHELL_EXPORT))
PROGRAMSOBJ1 = $(shell make -sC $(PROGRAMSDIR) -f Make.obj getobj $(SHELL_EXPORT))
PROGRAMSOBJ = $(subst programs/src/main_halidi.o,,$(PROGRAMSOBJ1)) # remove unwanted object
TESTOBJ = $(shell make -sC $(TESTDIR) -f Make.obj getobj $(SHELL_EXPORT))

OBJECTS += $(COREOBJ) $(FILEOBJ) $(INTEGRATOROBJ) $(MODELOBJ)

DIRLIST = $(COREDIR) $(FILEDIR) $(INTEGRATORDIR) $(MODELDIR) $(PROGRAMSDIR)
INCLUDE_PATH = $(foreach dir,$(DIRLIST),-I$(dir)/inc)

#=======================================
# Compiler options

RELEASE_FLAGS = -DNDEBUG -DBOOST_DISABLE_ASSERTS -O3

CLANG_CXXFLAGS = -Weverything -Wno-c++98-compat -Wno-missing-prototypes
CLANG_CXXFLAGS += -stdlib=libc++
CLANG_LDFLAGS += -stdlib=libc++

GCC_CXXFLAGS += -Wall -Wextra -ansi -pedantic
GCC_CXXFLAGS += -Wshadow -Wfloat-equal -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual
GCC_CXXFLAGS += -Wmissing-include-dirs -Wconversion -Wwrite-strings
GCC_CXXFLAGS += -Wold-style-cast -Wpointer-arith -Wcast-qual -Wcast-align
GCC_CXXFLAGS += -Wswitch-enum -Wundef -Wredundant-decls -Wstrict-null-sentinel -Wwrite-strings
GCC_CXXFLAGS += -Wunreachable-code -Wno-unknown-pragmas
GCC_CXXFLAGS += -Weffc++
GCC_LDFLAGS = 

export CPLUS_INCLUDE_PATH=/opt/local/include

ifeq ($(IS_GCC),0)
  CXXFLAGS = $(GCC_CXXFLAGS)
  LDFLAGS += $(GCC_LDFLAGS)
else ifeq ($(IS_CLANG),0)
  CXXFLAGS = $(CLANG_CXXFLAGS)
  LDFLAGS += $(CLANG_LDFLAGS)
else
  $(error Unkown compiler. Got '$(CXX_VERSION)')
endif

# We are using C++11 !
CXXFLAGS += -std=c++11

# Add defines
CXXFLAGS += $(foreach d,$(DEFINES),-D$(d))

# Export useful variables
export CXX
export CXXFLAGS
export LDFLAGS

LDLIBS := $(LIBDIR)/$(shell $(MAKE) -s CXX=$(CXX) \
	IS_GCC=$(IS_GCC) \
	IS_CLANG=$(IS_CLANG)\
	 CXX_VERSION="$(CXX_VERSION)" \
	-C $(LIBDIR) get_name)
build_libs := $(shell if [[ -e "$(LDLIBS)" ]]; then echo 0;else echo 1;fi)

# ============================================================
# objects targets

all: build_clang

.PHONY: build_clang
build_clang:
	$(MAKE) -C clang-build

.PHONY: build_gcc
build_gcc:
	$(MAKE) -C gcc-build


.PHONY: core
core:
	$(MAKE) -C $(COREDIR) -f Make.obj

.PHONY: file
file:
	$(MAKE) -C $(FILEDIR) -f Make.obj

.PHONY: integrators
integrators:
	$(MAKE) -C $(INTEGRATORDIR) -f Make.obj

.PHONY: model
model:
	$(MAKE) -C $(MODELDIR) -f Make.obj

.PHONY: programs
programs:
	$(MAKE) -C $(PROGRAMSDIR) -f Make.obj

# ============================================================
# Halidi model

# .PHONY: halidi_debug
# halidi_debug: OBJECTS += programs/src/main_halidi.o
# halidi_debug: objects_debug
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(HALIDI_EXEC)
# 	$(FINISH_PRINT)

# .PHONY: halidi
# halidi: CXXFLAGS += $(RELEASE_FLAGS)
# halidi: OBJECTS += programs/src/main_halidi.o
# halidi: objects
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(HALIDI_EXEC)
# 	$(FINISH_PRINT)

# ============================================================
# Koenigsberger model

# .PHONY: koenigs_debug
# koenigs_debug: OBJECTS += programs/src/main_koenigsberger.o
# koenigs_debug: objects_debug
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(KOENIGS_EXEC)
# 	$(FINISH_PRINT)

# .PHONY: koenigs
# koenigs: CXXFLAGS += $(RELEASE_FLAGS)
# koenigs: OBJECTS += programs/src/main_koenigsberger.o
# koenigs: objects
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(KOENIGS_EXEC)
# 	$(FINISH_PRINT)

# ============================================================
# Bifurcation diagram

# .PHONY: biffdiag
# #biffdiag: CXXFLAGS += $(RELEASE_FLAGS)
# biffdiag: CXXFLAGS += -g
# biffdiag: OBJECTS += programs/src/main_bifuraction_diagram.o
# biffdiag: objects
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(BIFDIAG_EXEC)
# 	$(FINISH_PRINT)

# ============================================================
# Main program

.PHONY: main_debug
main_debug: OBJECTS += $(PROGRAMSOBJ)
main_debug: objects_debug boostlib
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(MAIN_EXEC)
	$(FINISH_PRINT)

.PHONY: main
main: CXXFLAGS += $(RELEASE_FLAGS)
main: OBJECTS += $(PROGRAMSOBJ)
main: objects boostlib
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(MAIN_EXEC)
	$(FINISH_PRINT)

# .PHONY: koenigs
# koenigs: CXXFLAGS += $(RELEASE_FLAGS)
# koenigs: OBJECTS += programs/src/main_koenigsberger.o
# koenigs: objects
# 	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $(KOENIGS_EXEC)
# 	$(FINISH_PRINT)

# ============================================================
# Test

.PHONY: test
test: CXXFLAGS += -DTEST
test: CXXFLAGS := $(subst -Wshadow,,$(CXXFLAGS))
test: 
	$(MAKE) -C test -f Make.obj

test_clean: 
	$(MAKE) -C test -f Make.obj clean

# ============================================================
# Objects

.PHONY: objects
objects: core file integrators model programs

.PHONY: objects_debug
objects_debug: CXXFLAGS += -g
objects_debug: objects

# ============================================================
# Test

.PHONY: usage
usage:
	@echo "Available targets :"
#	@echo -e "  - biffdiag        \tBuild bifuraction diagram program"
	@echo -e "  - boostlib        \tBuild Boost libraries required by program"
	@echo -e "  - boostlib_force  \tSame as above, but force rebuild"
	@echo -e "  - clean           \tClean all objects file"
	@echo -e "  - clean_doc       \tDocumentation generated files"
	@echo -e "  - doc             \tGenerate documentation using Doxygen"
#	@echo -e "  - halidi[_debug]  \tBuild simulation program using the Halidi model"
#	@echo -e "  - koenigs[_debug] \tBuild simulation program using the Koenigsberger model"
	@echo -e "  - main[_debug]    \tBuild main simulation program"
	@echo -e "  - mrproper        \tClean all objects file + all executable + documentation"
	@echo -e "  - objects         \tCompile all source files to objects files
	@echo -e "  - test            \tBuild test program"
	@echo "  	NB: #_debug targets compile in debug mode"

# ============================================================
# Boost libraries required by main program

boostlib:
ifeq (1,$(build_libs))
	$(MAKE) -C $(LIBDIR)
endif

boostlib_force:
	$(MAKE) -C $(LIBDIR)

# ============================================================
# Other

.PHONY: doc
doc:
	doxygen $(DOCDIR)/Doxyfile


%.o: %.cpp
	$(CXX) -c $< -o $@ $(INCLUDE_PATH) $(CXXFLAGS) 

# For Makefile-debug purposes
print:
	@echo "$(CXX)"
	@echo "$(CXXFLAGS)"
	@echo "$(INCLUDE_PATH)"
	@echo "$(OBJECTS) $(PROGRAMSOBJ)"
	@echo "$(LDLIBS)"
	@echo "$(LDFLAGS)"
	@echo "$(build_libs)"
	@echo "$(OBJDIR)"

# ============================================================

FINISH_PRINT= @echo "********** Compilation finished ! **********"

.PHONY: clean
clean:
	$(MAKE) -C $(COREDIR) -f Make.obj clean
	$(MAKE) -C $(INTEGRATORDIR) -f Make.obj clean
	$(MAKE) -C $(FILEDIR) -f Make.obj clean
	$(MAKE) -C $(MODELDIR) -f Make.obj clean
	$(MAKE) -C $(PROGRAMSDIR) -f Make.obj clean
	$(MAKE) -C $(TESTDIR) -f Make.obj clean

.PHONY: clean_deps
clean_deps:
	$(MAKE) -C $(COREDIR) -f Make.obj clean_deps
	$(MAKE) -C $(INTEGRATORDIR) -f Make.obj clean_deps
	$(MAKE) -C $(FILEDIR) -f Make.obj clean_deps
	$(MAKE) -C $(MODELDIR) -f Make.obj clean_deps
	$(MAKE) -C $(PROGRAMSDIR) -f Make.obj clean_deps
	$(MAKE) -C $(TESTDIR) -f Make.obj clean_deps

.PHONY: clean_doc
clean_doc:
	rm -rf $(DOCDIR)/html
	rm -rf $(DOCDIR)/latex

.PHONY: mrproper
mrproper : clean clean_deps
	rm -f $(MAIN_EXEC)
	rm -f $(HALIDI_EXEC)
	rm -f $(KOENIGS_EXEC)
	rm -f $(TEST_EXEC)

