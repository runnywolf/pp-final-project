all: main.out

main.out: main.cpp
	rm -f main.out
	g++ -O3 -mavx2 -pthread -o main.out main.cpp

clean:
	rm -f main.out
