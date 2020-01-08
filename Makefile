CC=gcc
MPICC=mpicc
LIB_ALGO_NUM_PATH=/home/cisd-faverge/algonum/lib
CFLAGS=-O0 -g --std=gnu11 -fopenmp -I./include/ #-Wl,-rpath,$(LIB_ALGO_NUM_PATH)
LDFLAGS=-lm #-L$(LIB_ALGO_NUM_PATH) -lalgonum
SRC=util.c ddot.c dgemm.c daxpy.c dscal.c dgemv.c dger.c dgetrf.c dtrsm.c dgemm_omp.c dgetrf-omp.c dgemm-tile.c pdgetrf.c
OBJ=$(patsubst %.c, %.o, $(SRC))
TST=driver.c

all:lib/libmyblas.a

test:lib/libmyblas.a driver.o
	${CC} ${CFLAGS} -o test ${OBJ} driver.o -L./lib/ -lmyblas ${LDFLAGS}
	./test

%-perf:%-perf.o lib/libmyblas.a
	${CC} ${CFLAGS} -o tst/$@ $< -L./lib/ -lmyblas ${LDFLAGS} #-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

#-DMKL_ILP64 -m64 -I${MKLROOT}/include
%-perf.o:tst/%-perf.c
	${CC} ${CFLAGS} -c $< ${LDFLAGS}

p%-perf:p%-perf.o lib/libmyblas.a
	${MPICC} ${CFLAGS} -o tst/$@ $< -L./lib/ -lmyblas ${LDFLAGS} -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

p%-perf.o:tst/p%-perf.c
	${MPICC} ${CFLAGS} -DMKL_ILP64 -m64 -I${MKLROOT}/include -c $< ${LDFLAGS}

run-%-test: %-test
	./tst/$^

tests: run-dgemm-test

%-test:%-test.o lib/libmyblas.a
	${MPICC} ${CFLAGS} -o tst/$@ $< -L./lib/ -lmyblas ${LDFLAGS} -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

#-DMKL_ILP64 -m64 -I${MKLROOT}/include
%-test.o:tst/%-test.c
	${CC} ${CFLAGS} -DMKL_ILP64 -m64 -I${MKLROOT}/include -c $< ${LDFLAGS}

%.o:src/%.c
	${CC} ${CFLAGS} -c $< ${LDFLAGS}

p%.o:src/p%.c
	${MPICC} ${CFLAGS} -c $< ${LDFLAGS}

lib/libmyblas.a:$(OBJ)
	mkdir -p lib/
	ar crs lib/libmyblas.a $^

clean:
	rm *.o lib/libmyblas.a
