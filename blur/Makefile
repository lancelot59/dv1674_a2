# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -pg -Wunused -Wall -Wunused -pthread

run: clean all test

all: clean blur

test: testclean target1 target2 target3 target4

test1: clean all testclean target4

callgrind:
	sleep 5
	valgrind --tool=callgrind ./blur 15 data/im3.ppm outputfiles/im3.ppm

target1:
	./blur 15 data/im1.ppm outputfiles/im1.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im1_gprof.txt
	sleep 1
	rm gmon.out

target2:
	./blur 15 data/im2.ppm outputfiles/im2.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im2_gprof.txt
	sleep 1
	rm gmon.out

target3:
	./blur 15 data/im3.ppm outputfiles/im3.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im3_gprof.txt
	sleep 1
	rm gmon.out

target4:
	./blur 15 data/im4.ppm outputfiles/im4.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im4_gprof.txt
	sleep 1
	rm gmon.out

memtest: clean blur
	rm outputfiles/*.ppm massif*
	valgrind --tool=massif --massif-out-file=massif_im1.out ./blur 15 data/im1.ppm outputfiles/im1.ppm
	sleep 10
	valgrind --tool=massif --massif-out-file=massif_im2.out ./blur 15 data/im2.ppm outputfiles/im2.ppm
	sleep 10
	valgrind --tool=massif --massif-out-file=massif_im3.out ./blur 15 data/im3.ppm outputfiles/im3.ppm
	sleep 10
	valgrind --tool=massif --massif-out-file=massif_im4.out ./blur 15 data/im4.ppm outputfiles/im4.ppm

threshold: matrix ppm filters threshold.cpp
	$(CXX) $(CXXFLAGS) threshold.cpp matrix.o ppm.o filters.o -o threshold

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

clean:
	rm -rf blur threshold *.ppm *.o *.dSYM 2> /dev/null

testclean:
	rm outputfiles/*
