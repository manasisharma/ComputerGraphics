all: graphics_hw3

graphics_hw3: main.o main.o
	g++ main.o main.o -o graphics_hw3 -lglut -lGL -lGLU -lm

main.o: main.cpp
	g++ -c main.cpp

main.o: main.cpp 
	g++ -c main.cpp

clean: 
	rm -rf *.o graphics_hw3
