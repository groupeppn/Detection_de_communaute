SRCS=coded.cpp
SRCP=codedp.cpp
CC=g++ 
CFLAGS= -O3  -pg  
OPENMP_FLAGS=-fopenmp

codes: $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o codes

codep: $(SRCP)
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) $(SRCP) -o codep



serial:
	make codes

parallel:
	make codep


clean:
	rm -f code?

all:
	make clean
	make serial
	make parallel

