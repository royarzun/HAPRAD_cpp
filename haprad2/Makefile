PROG     = haprad20
EXE      = $(PROG).exe
TEST     = test.exe
SOURCES  = $(wildcard *.[f])
INC      = $(wildcard *.[inc])
OBJ      = rcdat.o fhaprad.o ihaprad.o  haprad_utils.o exclusive_model.o h4.o h3.o pkhff.o init_pdf.o semi_inclusive_model.o
OBJTEST  = rcdat_test.o fhaprad_test.o ihaprad.o  haprad_utils.o exclusive_model.o h4.o h3.o pkhff.o init_pdf.o semi_inclusive_model.o


CERN_ROOT=/apps/cernlib/i386_fc8/2005/
CERNLIBS = -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib
LIBS     = -L$(CERN_ROOT)/lib $(CERNLIBS)

COMPILER = gfortran

F77OPT  = -march=athlon64 -DLinux -fno-automatic -ffixed-line-length-none \
-fdollar-ok -fno-second-underscore -Iinclude

%.o: %.f $(INC)
	$(COMPILER) $(F77OPT) -c $< -o $@

.PHONY:
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(COMPILER) $(F77OPT) -o $(EXE) $(OBJ) $(LIBS)

.PHONY:
test:   $(TEST) 
	
$(TEST): $(OBJTEST) $(INC)
	$(COMPILER) $(F77OPT) -o $(TEST) $(OBJTEST) $(LIBS)

.PHONY:
clean:
	rm -f  *.o core

.PHONY:
distclean:
	rm -f $(EXE) $(TEST) res.dat test.dat

