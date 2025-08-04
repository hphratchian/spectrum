#
# This is a simple makefile for building spin-squared calculation code.
#
MQCDir       = $(mqcinstall)
MQCMODS      = $(MQCDir)/GNU/mod
MQCLIB       = $(MQCDir)/GNU/lib
LIBS         = -llapack -lblas -L$(MQCLIB)
F03Flags     = 
RunF         = gfortran -std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp
#RunF         = pgfortran -i8 -r8 -Mallocatable=03
#
#
# The 'all' rule.
#
all: spectrum.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
#hphmqc_integrals.mod: mqc_integrals.F03 $(MQCLIB)/libmqc.a
#hph	$(RunF) -I$(MQCMODS) -c $*.F03 $(MQCLIB)/libmqc.a

%.o: %.f90
	$(RunF) -J$(MQCMODS) -c $*.f90

#%.o: %.f03
#	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.f03

#
# Generic rule for building general executable program (*.exe) from a standard
# f90 source (*.f90) file.
#

%.exe: %.f03 spectrum_mod.f03 $(MQCLIB)/libmqc.a
	$(RunF) $(LIBS) $(Prof) -J$(MQCMODS) -o $*.exe $*.f03 $(MQCLIB)/libmqc.a

#
# Some clean rules.
#
cleanGTests:
	(cd GTests ; rm -f *.chk ; rm -f *.mat ; rm -f *.log)	
