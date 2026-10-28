CXX := c++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Iinclude

SRCS := src/bigint.cpp src/sha256.cpp src/ec.cpp src/ct.cpp src/ecdsa.cpp src/main.cpp
HDRS := include/bigint.hpp include/sha256.hpp include/ec.hpp include/ct.hpp include/ecdsa.hpp

.PHONY: all clean run

all: lab3

lab3: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

run: lab3
	./lab3

clean:
	rm -f lab3
