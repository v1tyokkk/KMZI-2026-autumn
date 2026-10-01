CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
SRCS      = gf8.cpp aes128.cpp gcm.cpp main.cpp

.PHONY: all clean

all: lab1

lab1: $(SRCS) gf8.hpp aes128.hpp gcm.hpp
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

clean:
	rm -f lab1
