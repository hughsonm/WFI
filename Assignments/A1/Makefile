# Usage:
# make # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all clean

CC = clang # compiler to use
CXX = clang++

CXXFLAGS = -O3 -march=native -msse4 -ffast-math -funroll-loops
LIBS = -lgmsh -lboost_math_C99
LDFLAGS = 

LINKERFLAG = -lm

SRCS := $(wildcard *.cpp)
BINS := $(SRCS:%.cpp=%)

# My A1 program depends on A1.o
A1: A1.o
	@echo "Building application A1"
	${CXX} ${CXXFLAGS} ${LIBS} ${LDFLAGS} -o A1.out A1.o #Build A1 from its dependencies

# My A1.o file depends on A1.cpp
A1.o: A1.cpp
	@echo "Building object file A1"
	${CXX} ${CXXFLAGS} -c -o A1.o  A1.cpp
# CXX and CXXFLAGS are special variables. make will use them to
# build A1.o. make is smart! I just need to tell it dependencies

clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BINS}
	rm -rvf *.out ${BINS}
