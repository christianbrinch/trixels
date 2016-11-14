shell = /bin/sh
FC = gfortran -m64 -I. -O3
CC = gcc -std=gnu99 -m64 -I. -O3
FL = gfortran -std=gnu99 -m64 -I. -O3

default: program

.SUFFIXES: .f90 .f .c .o
.c:
	$(CC) -c -o $*.o $<
.f.o:
	$(FC) -c -o $*.o $<

source:
	$(CC)  -unroll -c wrapper.c
	$(FC)  -unroll -c fourier.f

program: fourier.o wrapper.o
	$(FL) -o trixels.x *.o

clean:
	rm -f *.o *.x a.out
