# Usage:
# make # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all clean

CC = clang # compiler to use
CXX = clang++

CXXFLAGS = -O0
LIBS = -lgmsh
LDFLAGS = -L/usr/local/lib

LINKERFLAG = -lm

SRCS := $(wildcard *.cpp)
BINS := $(SRCS:%.cpp=%)

# My A1 program depends on A1.o
A1: A1.o
	${CXX} ${CXXFLAGS} ${LIBS} ${LDFLAGS} -o A1 A1.o #Build A1 from its dependencies

# My A1.o file depends on A1.cpp
A1.o: A1.cpp
	${CXX} ${CXXFLAGS} -c -o A1.o  A1.cpp
# CXX and CXXFLAGS are special variables. make will use them to
# build A1.o. make is smart! I just need to tell it dependencies


#edit:
#	${CXX} ${FLAGS} -o A1.out A1.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BINS}
	rm -rvf *.out ${BINS}
