CFLAGS = -Wall -g -Wno-vla 
CXX = clang++

all: build/main

build/main: main.cpp kenken.cpp solve.cpp
	mkdir -p build
	$(CXX) $(CFLAGS) $^ -o $@

clean:
	rm -rf build
