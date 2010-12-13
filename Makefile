SHELL = /bin/bash

.DELETE_ON_ERROR:

.SECONDARY: $(DICT_CLASS)
.PHONY: all lib clean distclean checkdirs

##############################################################################
# FORTRAN
##############################################################################

FEXE      := haprad20.exe
FTEST     := test.exe

INC      := $(wildcard include/*.inc)

FOBJ      := rcdat.o fhaprad.o ihaprad.o haprad_utils.o exclusive_model.o  \
	    h4.o h3.o pkhff.o init_pdf.o semi_inclusive_model.o

OBJ_FTEST  := rcdat_test.o fhaprad_test.o ihaprad.o haprad_utils.o          \
	    exclusive_model.o h4.o h3.o pkhff.o init_pdf.o                \
	    semi_inclusive_model.o

CERNLIBS  := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib
FLIBS      := $(CERNLIBS)

FCC  := gfortran

F77OPT    := -march=i686 -DLinux -fno-automatic -ffixed-line-length-none \
             -fdollar-ok -fno-second-underscore -Iinclude

##############################################################################
# C++
##############################################################################

ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTCINT    := rootcint

ifeq "$(DEBUG)" ""
CXXFLAGS += -O2
else
CXXFLAGS += -g
endif

ifeq ($(findstring Linux,$(OS_NAME)),Linux)
CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
endif

LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)
SOFLAGS   = -Wl,-soname,$(notdir $@) -shared

INCLUDES  = -I.
LIBS      =


SLIB_DIR := slib
OBJ_DIR  := .obj
DICT_DIR := .dict
DEP_DIR  := .dep

SRC_CLASS  := $(wildcard *.cxx)
SRC_DEP    := $(addprefix $(DEP_DIR)/,$(SRC_CLASS:.cxx=.d))
SRC_FOBJ    := $(addprefix $(OBJ_DIR)/,$(SRC_CLASS:.cxx=.o))
DICT_CLASS := $(addprefix $(DICT_DIR)/,$(SRC_CLASS:.cxx=Dict.cxx))
DICT_FOBJ   := $(addprefix $(OBJ_DIR)/,$(SRC_CLASS:.cxx=Dict.o))

SH_LIB     := libTRadCor.so

##############################################################################
# Rules
##############################################################################

all: $(FEXE) lib

##############################################################################

lib: checkdirs $(SLIB_DIR)/$(SH_LIB)

include Makefile_depends


$(SLIB_DIR)/$(SH_LIB): $(SRC_FOBJ) $(DICT_FOBJ)
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@


$(DICT_DIR)/%Dict.cxx: %.h
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(DICT_DIR) || mkdir -p $(DICT_DIR)
endif
	$(ROOTCINT) -f $@ -c $(<F) $(<F:%.h=%LinkDef.h)


$(OBJ_DIR)/%Dict.o: $(DICT_DIR)/%Dict.cxx
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)
endif
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


$(OBJ_DIR)/%.o: %.cxx
ifneq (,$(filter %.o,$(MAKECMDGOALS)))
	@test -d $(DEP_DIR) || mkdir -p $(DEP_DIR)
	@test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)
endif
	@$(call make-depend,$<,$@,$(@F:.o=.d))
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


checkdirs: $(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR)

$(SLIB_DIR) $(OBJ_DIR) $(DICT_DIR) $(DEP_DIR):
	@mkdir -p $@

##############################################################################

$(FEXE): $(FOBJ) $(INC)
	$(FCC) $(F77OPT) -o $(FEXE) $(FOBJ) $(FLIBS)

test: $(FTEST)

$(FTEST): $(OBJ_FTEST) $(INC)
	$(FCC) $(F77OPT) -o $(FTEST) $(OBJ_FTEST) $(FLIBS)


%.o: %.f $(INC)
	$(FCC) $(F77OPT) -c $< -o $@

##############################################################################

clean:
	@rm -f  *.o core
	@rm -rf $(OBJ_DIR) $(DICT_DIR)

distclean: clean
	@rm -f $(FEXE) $(FTEST) res.dat test.dat
	@rm -rf $(SLIB_DIR) $(DEP_DIR)
