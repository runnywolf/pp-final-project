all: main.out

main.out: main.cpp
	rm -f main.out
	g++ -O3 -g -mavx2 -pthread -o main.out main.cpp

clean:
	rm -f main.out

record:
	perf record -o perf.data ./main.out
