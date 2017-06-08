######################################################################
## This is the make file to compile LASER and TRACE programs
## Assuming the LAPACK and BLAS libraries were installed in /usr/lib64
## Assuming the GSL and Armadillo libraries were installed in ~/usr/lib
## Statically link to libgsl.a, libgslcblas.a, liblapack.a, libblas.a
##
## Type "make laser" to compile LASER
## Type "make trace" to compile TRACE
## Type "make" to compile both LASER and TRACE
#######################################################################

all: laser trace

LOCAL = /home/chaolong/usr
## Change $(LOCAL) if GSL and Armadillo were installed in a different directory

G++FLAG = -O3 -fopenmp

laser: laser.o
	g++ -static -o laser laser.o \
	/usr/lib/libgsl.a \
	/usr/lib/libgslcblas.a \
	/usr/lib/liblapack.a \
	$(LOCAL)/lib/libopenblas.a \
	-lgfortran -lpthread

laser.o: laser.cpp
	g++ $(G++FLAG) -c laser.cpp -I /usr/include -I $(LOCAL)/include  
	

trace: trace.o
	g++ -static -o trace trace.o \
	/usr/lib/libgsl.a \
	/usr/lib/libgslcblas.a \
	/usr/lib/liblapack.a \
	$(LOCAL)/lib/libopenblas.a \
	-lgfortran -lpthread

trace.o: trace.cpp
	g++ $(G++FLAG) -c trace.cpp -I /usr/include -I $(LOCAL)/include  


clean:
	rm *.o 

