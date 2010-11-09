EXE      := haprad20.exe
TEST     := test.exe

INC      := $(wildcard include/*.inc)

OBJ      := rcdat.o fhaprad.o ihaprad.o haprad_utils.o exclusive_model.o  \
	    h4.o h3.o pkhff.o init_pdf.o semi_inclusive_model.o

OBJTEST  := rcdat_test.o fhaprad_test.o ihaprad.o haprad_utils.o          \
	    exclusive_model.o h4.o h3.o pkhff.o init_pdf.o                \
	    semi_inclusive_model.o


CERNLIBS  := -lpdflib804 -lmathlib -lphtools -lpacklib -lkernlib -lpawlib
LIBS      := $(CERNLIBS)

COMPILER  := gfortran

F77OPT    := -march=i686 -DLinux -fno-automatic -ffixed-line-length-none \
             -fdollar-ok -fno-second-underscore -Iinclude

all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(COMPILER) $(F77OPT) -o $(EXE) $(OBJ) $(LIBS)

test: $(TEST)
	
$(TEST): $(OBJTEST) $(INC)
	$(COMPILER) $(F77OPT) -o $(TEST) $(OBJTEST) $(LIBS)


%.o: %.f $(INC)
	$(COMPILER) $(F77OPT) -c $< -o $@


clean:
	rm -f  *.o core

distclean: clean delete

delete:
	rm -f $(EXE) $(TEST) res.dat test.dat
