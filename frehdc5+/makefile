# makefile for the implicit advdiff solver
HOME=laspack
#include $(HOME)

CC = mpicc

all:
	$(CC) bathymetry.c configuration.c fileio.c initialize.c map.c \
		mpifunctions.c nsfunctions.c nssolve.c scalar.c subgrid.c utilities.c \
	  $(HOME)/eigenval.c $(HOME)/errhandl.c $(HOME)/factor.c \
		$(HOME)/itersolv.c $(HOME)/matrix.c $(HOME)/mlsolv.c \
		$(HOME)/operats.c $(HOME)/precond.c $(HOME)/qmatrix.c $(HOME)/vector.c \
		$(HOME)/rtc.c FREHD.c -lm -o runthis.o

	#$(CC) bathymetry.c configuration.c fileio.c map.c mpifunctions.c utilities.c FREHD.c -o runthis.o
