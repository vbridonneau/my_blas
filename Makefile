CC=gcc
CFLAGS=-std=c99
LDFLAGS=-lm
SRC=util.c ddot.c dgemm.c daxpy.c dscal.c dgemv.c
OBJ=$(patsubst %.c, %.o, $(SRC))
TST=driver.c

all:libmyblas.a

test:libmyblas.a driver.o
	${CC} ${CFLAGS} -o test ${OBJ} driver.o -L. -lmyblas ${LDFLAGS}
	./test

%-perf:libmyblas.a %-perf.o
	${CC} ${CFLAGS} -o $@ $^ -L. -lmyblas ${LDFLAGS}
	sh ./$@.sh

%.o:%.c
	${CC} ${CFLAGS} -c $< ${LDFLAGS}

libmyblas.a:$(OBJ)
	ar crs libmyblas.a $^

clean:
	rm *.o libmyblas.a
