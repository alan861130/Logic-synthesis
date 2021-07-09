.PHONY: all clean

all:	map

decompose:	decompose.cpp
	$(CXX) decompose.cpp -o decompose 

map:	map.cpp
	$(CXX) map.cpp -o map 

clean:
	$(RM) map decompose output.blif
