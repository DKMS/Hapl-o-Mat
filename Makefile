CXX = g++
#CXXFLAGS = -W -Wall -g -Wextra -pedantic -std=c++14 -I include
CXXFLAGS = -Wall -march=native -Ofast -std=c++14 -I include

OBJ = $(patsubst %.cc, %.o, $(wildcard src/*.cc))

haplomat: $(OBJ)
	$(CXX) ${CXXFLAGS} $(OBJ) -o $@ 

%.o: %.cc %.h
	$(CXX) ${CXXFLAGS} -c $< -o $@

clean:
	$(RM) -f src/*.o
