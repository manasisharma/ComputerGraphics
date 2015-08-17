CC=g++
FLAGS=-lGL -lglut
SRC=p2.cpp p5.cpp part1_final.cpp
OBJ=$(SRC:.cpp=.o)
all: $(SRC) $(OBJ)
.cpp.o:
	$(CC) $< $(FLAGS) -o $@
clean:
	rm $(OBJ)
