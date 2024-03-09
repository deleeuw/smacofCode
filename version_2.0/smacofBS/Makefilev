CFLAGS = -O2 -Wall -Wextra -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -pedantic -std=c99

sources = \
	smacofUtils.c \
	smacofMatrixMake.c \
	smacofSetup.c \
	smacofLinearAlgebra.c \
	smacofNorm.c \
	smacofInitial.c \
	smacofBasis.c \
	smacofCCD.c \
	smacofPrint.c \
	smacofParameters.c \
	smacofCore.c 

objects = \
	smacofUtils.o \
	smacofMatrixMake.o \
	smacofSetup.o \
	smacofLinearAlgebra.o \
	smacofNorm.o \
	smacofInitial.o \
	smacofBasis.o \
	smacofCCD.o \
	smacofPrint.o \
	smacofParameters.o \
	smacofCore.o 
	
drivers = smacofBS.c \
	smacofBSEngine.c

includes = smacof.h

all: smacofBS

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofBS.a $(objects)
	mv -f libsmacofBS.a /usr/local/lib
	rm -f *.o

smacofBS: clib
	$(CC) -o smacofBS $(CFLAGS) $(drivers) -lsmacofBS
	rm -f /usr/local/bin/smacofBS
	mv -f smacofBS /usr/local/bin
	
smacofBSDebug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofBSDebug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofBS.so smacofBSEngineR.c $(sources)

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofBS smacofBSDebug \
		libsmacofBS.a libsmacofBS.so libsmacofBS.0.dylib libsmacofBS.dylib

