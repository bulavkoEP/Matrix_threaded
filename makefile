CFLAGS=-fsanitize=address -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wmissing-declarations -lm -O3

LDFLAGS=-fsanitize=address

all: prog

prog: main.o io.o functions.o Cholesky.o
	g++ main.cpp io.o functions.o Cholesky.o -o main $(LDFLAGS)

io.o: io.cpp io.h
	g++ -c $(CFLAGS) io.cpp

functions.o: functions.cpp functions.h
	g++ -c $(CFLAGS) functions.cpp

Gauss.o: Cholesky.cpp Cholesky.h
	g++ -c $(CFLAGS) Cholesky.cpp

clean:
	rm -rf *.o prog