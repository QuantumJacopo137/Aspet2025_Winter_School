F90=gfortran
##LDFLAGS=-L/usr/lib64 -llapack

all:
	$(F90) -c Heis.f Det.f Factorial.f ReadInput.f Integ_to_Bit.f ValCoupl.f Diasym.f Identify_Ms.f QEqual.f
	$(F90) -o Heis Heis.o Det.o Factorial.o ReadInput.o Integ_to_Bit.o ValCoupl.o Diasym.o Identify_Ms.o QEqual.o /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#	$(F90) -o Heis Heis.o Det.o Factorial.o ReadInput.o Integ_bool.o ValCoupl.o Diasym.o Identify_Ms.o QEqual.o /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3 /usr/lib/x86_64-linux-gnu/blas/libblas.so.3
clean:
	rm -rf *.o
	rm -f Deter
