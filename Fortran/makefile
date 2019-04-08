
FemPear:inout2.o matrix_solve.o minpackfunctions.o hybrd.o FemPear.o
	$ gfortran -Wall -o FemPear inout2.o  matrix_solve.o minpackfunctions.o hybrd.o FemPear.o
	./FemPear > result.out -Wall -Wtabs


inout2.o:./inout2.f90
	$ gfortran -Wall -c ./inout2.f90

matrix_solve.o:./matrix_solve.f90
	$ gfortran -Wall -c ./matrix_solve.f90

minpackfunctions.o:./minpackfunctions.f90
	$ gfortran -Wall -c ./minpackfunctions.f90

hybrd.o:./hybrd.f90 minpackfunctions.o
	$ gfortran -Wall -c ./hybrd.f90


FemPear.o:FemPear.f90 inout2.o matrix_solve.o minpackfunctions.o hybrd.o
	$ gfortran -Wall -c ./FemPear.f90


clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M  *.d V*.inc *.vo \
	V*.f *.dbg album F.err
