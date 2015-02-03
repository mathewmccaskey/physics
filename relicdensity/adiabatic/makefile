fort = gfortran

flags = -O -ffixed-line-length-132

prog = adiabatic

.f.o:
	$(fort) -c $(flags) $*.f

obj = entropy.o romberg.o

all: clean jazz

clean:
	-rm *.o
	-rm $(prog).x

jazz: $(obj)
	$(fort) $(flags) $(obj) $(libs) -o $(prog).x