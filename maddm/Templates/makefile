all:
	-cd Source/DHELAS/ && make
	-cd Source/MODEL/ && make
	-cd matrix_elements/ && make
	-mv matrix_elements/libmatrix_elements.a lib/
	-ranlib lib/libmatrix_elements.a
	-cd src/ && make
	-mv src/maddm.x .

clean:
	-cd matrix_elements/ && make clean
	-cd src/ && make clean
