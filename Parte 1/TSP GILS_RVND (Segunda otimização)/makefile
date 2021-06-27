src = $(wildcard src/*.cpp)
obj = $(src:.cpp=.o)

CXXFLAGS = -std=c++1z -Wall -g -O2

tsp: $(obj)
	$(CXX) -g -o $@ $^ $(CXXFLAGS)
clean:

	rm -rf *o tsp
	rm -rf src/*.o
