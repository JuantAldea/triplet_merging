TARGET = thread
FLAGS = -std=c++11 -g  -Wall -Wextra -O0 -lpthread -lstdc++ -lm

all: gcc

clang:
	clang -o $(TARGET) thread.cpp $(FLAGS)

gcc: 
	g++ -o $(TARGET) thread.cpp $(FLAGS) 

clean:
	rm -rf thread

.PHONY : all
