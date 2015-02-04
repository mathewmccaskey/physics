prog = maddm
fort= gfortran

libdir = ../lib/
incdir = ../include

flags = -O -ffixed-line-length-132 -L$(libdir) -I$(incdir)

libs = -lmatrix_elements -ldhelas -lmodel

objs = init.o integrate.o interpolate.o maddm.o phasespace.o relic_canon.o relic_coupled.o tests.o

.f.o:
	$(fort) -c $(flags) $*.f

all: clean maddm

clean:
	-rm -f *.o
	-rm -f ../$(prog).x

maddm: $(objs)
	$(fort) $(flags) $(objs) $(libs) -o $(prog).x
