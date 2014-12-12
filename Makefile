CXX = g++
INC1 = include
INCDIRS = -I${INC1}
CXXFLAGS = -W -Wall -g -Wextra -pedantic -std=c++11 ${INCDIRS}
#CXXFLAGS = -Wall -pg -march=native -O2 -fno-rtti -std=c++11 ${INCDIRS}
#CXXFLAGS = -Wall -fopenmp -march=native -Ofast -fno-rtti -std=c++11 ${INCDIRS}

sources = $(wildcard *.cc)
OBJ = $(patsubst %.cc, %.o, $(wildcard *.cc))

main: $(OBJ)
	$(CXX) ${CXXFLAGS} -o $@ $(OBJ) 

%.d: %.cc %.h
	@ $(CXX) -MM $(CXXFLAGS) $< -o $@

%.o: %.cc %.h
	$(CXX) ${CXXFLAGS} -c $< -o $@

clean:
	$(RM) -f *.o *.d
-include $(sources:.cc=.d)

.PHONY: clean all
