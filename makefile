all: main.cpp tools.h fft.h
	g++ -std=c++11  main.cpp -g -o poly-mult