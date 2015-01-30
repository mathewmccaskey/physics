all: cleanlibs makelibs

cleanlibs:
	-rm *.a

makelibs: matlib

matlib:
	-cd src/ && make
	-mv src/libmat.a .
	-ranlib libmat.a
