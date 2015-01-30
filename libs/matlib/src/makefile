fort = ifort
home = $(shell pwd)
libmake = $(home)
libdir = $(libmake)
incdir = $(home)/include

flags = -g -O3 -fpp -axP -parallel -132 -ip -L$(libidr) -D$(shell uname -s)

.f.o:
	$(fort) -c $(flags) -I$(incdir) $*.f

obj = amoeba.o analyzeevent.o besselk.o bessi0.o bessi1.o bessk0.o bessk1.o bessk.o brent.o collidertools.o decaymother.o eigsrt.o hist.o interpolate.o jacobi.o maths.o matrixmath.o mt2.o odeint.o phasespace.o polint.o qromb.o qsimp.o random.o readevents.o romberg.o setgamma.o simpson.o spline.o splint.o trapzd.o vegas.o

libname = mat

all: cleanlib makematlib

cleanlib:
	-rm *.o

makematlib: $(obj)
	-ar r lib$(libname).a $(obj)
