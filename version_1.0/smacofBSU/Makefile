CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -pedantic -std=c99

sources = \
	smacofPrintReadU.c \
	smacofUtilsU.c \
  smacofMatrixUtils.c \
	smacofLinearAlgebra.c \
	smacofNormU.c \
	smacofInitial.c \
	smacofBases.c \
	smacofCCDU.c \
	smacofDerivativesU.c \
	smacofCoreU.c 

objects = \
	smacofPrintReadU.o \
	smacofUtilsU.o \
	smacofMatrixUtils.o \
	smacofLinearAlgebra.o \
	smacofNormU.o \
	smacofInitial.o \
	smacofBases.o \
	smacofCCDU.o \
	smacofDerivativesU.o \
	smacofCoreU.o 
	
drivers = smacofBSU.c \
	smacofBSUEngine.c

includes = smacofU.h

all: smacofBSU

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofBSU.a $(objects)
	mv -f libsmacofBSU.a /usr/local/lib
	rm -f *.o

smacofBSU: clib
	$(CC) -o smacofBSU $(CFLAGS) $(drivers) -lsmacofBSU
	rm -f /usr/local/bin/smacofBSU
	mv -f smacofBSU /usr/local/bin
	
smacofBSUDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofBSUDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofBSU.so smacofBSUEngineR.c $(sources)

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofBSU smacofBSUDebug \
		libsmacofBSU.a libsmacofBSU.so libsmacofBSU.0.dylib libsmacofBSU.dylib

