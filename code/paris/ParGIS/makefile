






../lib/libParGIS.a: ParGIS.o
	ar -crv ../lib/libParGIS.a ParGIS.o
ParGIS.o: ParGIS.cu ParGIS.h
	/usr/local/cuda/bin/nvcc -c ParGIS.cu -o ParGIS.o -lcudart -lcuda



#hello.a:
#	ar -crv libhello.a libhello.o

#hello.o: hello.cu hello.h
#	g++ -c hello.cpp -o libhello.o

