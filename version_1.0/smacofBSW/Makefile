CFLAGS = -O2 -Wall -Wextra -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -pedantic -std=c99

sources = \
	smacofPrintReadW.c \
	smacofUtilsW.c \
  	smacofMatrixUtils.c \
	smacofLinearAlgebra.c \
	smacofNormW.c \
	smacofInitial.c \
	smacofBases.c \
	smacofCCDW.c \
	smacofDerivativesW.c \
	smacofCoreW.c 

objects = \
	smacofPrintReadW.o \
	smacofUtilsW.o \
	smacofMatrixUtils.o \
	smacofLinearAlgebra.o \
	smacofNormW.o \
	smacofInitial.o \
	smacofBases.o \
	smacofCCDW.o \
	smacofDerivativesW.o \
	smacofCoreW.o 
	
drivers = smacofBSW.c \
	smacofBSWEngine.c

includes = smacofW.h

all: smacofBSW

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofBSW.a $(objects)
	mv -f libsmacofBSW.a /usr/local/lib
	rm -f *.o

smacofBSW: clib
	$(CC) -o smacofBSW $(CFLAGS) $(drivers) -lsmacofBSW
	rm -f /usr/local/bin/smacofBSW
	mv -f smacofBSW /usr/local/bin
	
smacofBSWDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofBSWDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofBSW.so smacofBSWEngine.c $(sources)

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofBSW smacofBSWDebug \
		libsmacofBSW.a libsmacofBSW.so libsmacofBSW.0.dylib libsmacofBSW.dylib

