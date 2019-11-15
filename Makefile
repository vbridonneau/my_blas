CC=gcc
LIB_ALGO_NUM_PATH=/home/cisd-faverge/algonum/lib
CFLAGS=-std=c99 -I./include/ -Wl,-rpath,$(LIB_ALGO_NUM_PATH)
LDFLAGS=-lm -L$(LIB_ALGO_NUM_PATH) -lalgonum
SRC=util.c ddot.c dgemm.c daxpy.c dscal.c dgemv.c dger.c
OBJ=$(patsubst %.c, %.o, $(SRC))
TST=driver.c

all:lib/libmyblas.a

test:lib/libmyblas.a driver.o
	${CC} ${CFLAGS} -o test ${OBJ} driver.o -L./lib/ -lmyblas ${LDFLAGS}
	./test

%-perf:lib/libmyblas.a %-perf.o
	${CC} ${CFLAGS} -o tst/$@ $^ -L./lib/ -lmyblas ${LDFLAGS}
	./tst/$@

%-perf.o:tst/%-perf.c
	${CC} ${CFLAGS} -c $< ${LDFLAGS}

%.o:src/%.c
	${CC} ${CFLAGS} -c $< ${LDFLAGS}

lib/libmyblas.a:$(OBJ)
	ar crs lib/libmyblas.a $^

clean:
	rm *.o lib/libmyblas.a
