CXX := g++
SRC := $(wildcard src/*.cpp) $(wildcard src/arboretum/*.cpp) 
OBJ := $(SRC:.cpp=.o)

CXXFLAGS := -Iinclude `pkg-config --cflags opencv` -O3 -g -fPIC -Wall -Wextra -std=c++11 -pedantic
LDFLAGS := -llmdb -lboost_system -lboost_filesystem  `pkg-config --libs opencv` -O3 -g

all: Halite libhalite.so

Halite: demo/Halite.o libhalite.a
	$(CXX) demo/Halite.o libhalite.a $(LDFLAGS) -o Halite

libhalite.so: $(OBJ)
	$(CXX) $(OBJ) -shared -o libhalite.so

libhalite.a: $(OBJ)
	ar -rcs libhalite.a $(OBJ)
demo:   Halite
	./Halite 1 1e-10 4 1 databases/12d.dat result12d.dat
clean:
	rm -f $(OBJ) demo/Halite.o

spotless: clean
	rm -f Halite

obj/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
