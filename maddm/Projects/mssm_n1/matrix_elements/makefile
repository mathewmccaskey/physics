fort= gfortran

libname = matrix_elements
incdir = ../include

flags = -O -ffixed-line-length-132 -I$(incdir)

objs = smatrix.o matrix_n1n1_uux.o matrix_n1n1_ddx.o matrix_n1n1_vevex.o matrix_n1n1_vtvtx.o matrix_n1n1_emep.o matrix_n1n1_ttx.o matrix_n1n1_bbx.o matrix_n1n1_zz.o matrix_n1n1_zh1.o matrix_n1n1_wpwm.o matrix_n1n1_h1h1.o matrix_n1n1_tamtap.o

.f.o:
	$(fort) -c $(flags) $*.f

all: make_lib

clean:
	-rm -f *.o
	-rm -f ../lib/lib$(libname).a

make_lib: $(objs)
	-ar r lib$(libname).a $(objs)
