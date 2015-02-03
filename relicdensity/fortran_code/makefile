fort = ifort

libdir = /Users/mccaskey/Research/libs/matlib/
incdir = $(libdir)/src/include

flags = -g -O3 -fpp -axP -parallel -ip -132 -L$(libdir) -I$(incdir)

libs = -lmat

prog = relic

.f.o:
	$(fort) -c $(flags) $*.f

#obj = relicdensity.f_OLD.o
obj = relicdensity.o

all: clean jazz

clean:
	-rm *.o
	-rm $(prog).x

jazz: $(obj)
	$(fort) $(flags) $(obj) $(libs) -o $(prog).x
