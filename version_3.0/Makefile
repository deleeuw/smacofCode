CFLAGS = -O2 -Wall -Wextra -Werror -Wnull-dereference -pedantic -std=c99
DFLAGS = -g3 -Wall -Wextra -Werror -Wnull-dereference -pedantic -std=c99

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

all: smacofBS3

ofiles: $(includes) $(sources)
	$(CC) -c $(CFLAGS) $(sources)

clib: ofiles 
	ar rcs libsmacofBS3.a $(objects)
	mv -f libsmacofBS3.a /usr/local/lib
	rm -f *.o

smacofBS3: clib
	$(CC) -o smacofBS3 $(CFLAGS) $(drivers) -lsmacofBS3
	rm -f /usr/local/bin/smacofBS3
	mv -f smacofBS3 /usr/local/bin
	
smacofBS3Debug: $(includes) $(sources) $(drivers)
	$(CC) -o smacofBS3Debug $(DFLAGS) $(drivers) $(sources)

rlib: ofiles
	R CMD SHLIB -o libsmacofBS3.so smacofBS3EngineR.c $(sources)

format: $(includes) $(sources) $(drivers)
	clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4}" \
	$(sources) $(includes) $(drivers)
	
clean:
	rm -f $(objects) ../*/*.o

pristine: clean
	rm -f smacofBS3 smacofBSDebug3 \
		libsmacofBS3.a libsmacofBS3.so libsmacofBS3.0.dylib libsmacofBS3.dylib

