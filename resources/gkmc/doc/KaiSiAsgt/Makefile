#SHELL=/bin/tcsh
CFLAGS = -Wall -g
INCLUDES = -I. -I/usr/include -I/afs/ipp/home/m/mow/soft/gsl-1.2.90/include
LIBS = -L/afs/ipp/home/m/mow/soft/gsl-1.2.90/lib -lgsl -lgslcblas -lm

CC= gcc

#Rule to create .o files from .c files
.c.o:
	 $(CC) -c $(CFLAGS) $(INCLUDES) $*.c

SiAsgt : SiAsgt.o allocarrays.o readinputs.o printinputs.o initmaterial.o \
	initparams.o bklalgorithm.o
	$(CC) SiAsgt.o  allocarrays.o readinputs.o printinputs.o \
	initmaterial.o initparams.o bklalgorithm.o $(LIBS) -o SiAsgt

SiAsgt.o: SiAsgt.c allocarrays.c readinputs.c printinputs.c \
	initparams.c bklalgorithm.c
allocarrays.o: allocarrays.c
readinputs.o: readinputs.c
printinputs.o: printinputs.c
bklalgorithm.o: bklalgorithm.c initmaterial.c
initmaterial.o: initmaterial.c
initparams.o: initparams.c

clean : 
	/bin/rm *.o *.xyz *.log SiAsgt

run :
	./SiAsgt
