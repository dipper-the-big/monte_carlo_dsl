#SHELL=/bin/tcsh
CFLAGS = -Wall -g -O2
INCLUDES = -I. -I/usr/include -I/usr/local/include
LIBS = -L/usr/local/lib -lgsl -lgslcblas -lm

CC= gcc

#Rule to create .o files from .c files
.c.o:
	 $(CC) -c $(CFLAGS) $(INCLUDES) $*.c


RDAlenya : RDAlenya.o readinputs.o printinputs.o initparams.o \
      kmcdriver.o allocstatvarmem.o initparticles.o posout.o \
      bklalgorithm.o applybc.o findcellnum.o makeajump.o \
      jumponsurface.o initgeometry.o reallocmem.o checkrecomb.o findinteract.o \
      interactfunc.o getdiff.o desorb.o postprocess.o diffcoeffeval.o statpostproc.o \
      dostats.o
	$(CC) RDAlenya.o readinputs.o printinputs.o initparams.o  \
      kmcdriver.o allocstatvarmem.o initparticles.o posout.o \
      bklalgorithm.o applybc.o findcellnum.o makeajump.o \
      jumponsurface.o initgeometry.o reallocmem.o checkrecomb.o findinteract.o \
      interactfunc.o getdiff.o desorb.o postprocess.o  diffcoeffeval.o statpostproc.o \
      dostats.o $(LIBS) -o RDAlenya

RDalenya.o: RDAleny.c readinputs.c printinputs.c initparams.c initgeometry.o \
       kmcdriver.c GlobalVars.h
readinputs.o: readinputs.c GlobalVars.h 
printinputs.o: printinputs.c GlobalVars.h
initparams.o: initparams.c GlobalVars.h GlobalGSL.h
initgeometry.o: initgeometry.c GlobalVars.h
kmcdriver.o: kmcdriver.c allocstatvarmem.c initparticles.c \
       bklalgorithm.c reallocmem.c posout.c initgeometry.c \
       postprocess.c statpostproc.c GlobalVars.h 
allocstatvarmem.o: allocstatvarmem.c GlobalVars.h 
initparticles.o: initparticles.c GlobalVars.h GlobalGSL.h
bklalgorithm.o: bklalgorithm.c \
       jumponsurface.c checkrecomb.c GlobalVars.h \
       GlobalGSL.h
applybc.o: applybc.c GlobalVars.h 
findcellnum.o: findcellnum.c GlobalVars.h 
makeajump.o: makeajump.c GlobalVars.h 
jumponsurface.o: jumponsurface.c applybc.c findcellnum.c makeajump.c \
       GlobalVars.h GlobalGSL.h
checkrecomb.o: checkrecomb.c findinteract.c interactfunc.c GlobalVars.h \
    GlobalGSL.h
findinteract.o: findinteract.c GlobalVars.h 
interactfunc.o: interactfunc.c getdiff.c GlobalVars.h 
getdiff.o: getdiff.c GlobalVars.h 
desorb.o: desorb.c GlobalVars.h
reallocmem.o: reallocmem.c GlobalVars.h 
postprocess.o: postprocess.c diffcoeffeval.c GlobalVars.h 
diffcoeffeval.o: diffcoeffeval.c getdiff.c GlobalVars.h 
statpostproc.o: statpostproc.c dostats.c GlobalVars.h GlobalGSL.h
dostats.o: dostats.c GlobalGSL.h
posout.o: posout.c GlobalVars.h

clean : 
	/bin/rm *.o RDAlenya *~

cleanout : 
	/bin/rm ../out/*.out ../out/*.log ../out/*.dat ../out/*.xyz

run :
	time ./RDAlenya
