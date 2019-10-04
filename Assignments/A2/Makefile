# Usage:
# make # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all clean

CC = clang # compiler to use
CXX = clang++

CXXFLAGS = -O3 -march=native -msse4 -ffast-math -funroll-loops -Wall
LIBS = -lgmsh -lboost_math_C99
LDFLAGS =

LINKERFLAG = -lm

SRCS := $(wildcard *.cpp)
BINS := $(SRCS:%.cpp=%)

# My A1 program depends on A1.o
A2: A2.o ../A1/mom_driver.o
	@echo "Building application A1"
	${CXX} ${CXXFLAGS} ${LIBS} ${LDFLAGS} -o A2.out A2.o ../A1/mom_driver.o #Build A1 from its dependencies

# My A1.o file depends on A1.cpp
A2.o: A2.cpp
	@echo "Building object file A1"
	${CXX} ${CXXFLAGS} -c -o A2.o A2.cpp
# CXX and CXXFLAGS are special variables. make will use them to
# build A1.o. make is smart! I just need to tell it dependencies

../A1/mom_driver.o: ../A1/mom_driver.cpp ../A1/mom_driver.h
	@echo "Building object file mom_driver"
	${CXX} ${CXXFLAGS} -c -o ../A1/mom_driver.o ../A1/mom_driver.cpp

clean:
	@echo "Cleaning up..."
	rm -rvf *.o ${BINS}
	rm -rvf *.out ${BINS}
