# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -pg -Wunused -Wall -Wunused -pthread

run: clean all test

all: clean blur blur_par

test:  target1 target2 target3 target4

target1:
	./blur_par 15 data/im1.ppm outputfiles_par/im1.ppm 8
	sleep 1
	gprof blur_par gmon.out > outputfiles_par/im1_gprof.txt
	sleep 1
	rm gmon.out

target2:
	./blur_par 15 data/im2.ppm outputfiles_par/im2.ppm 8
	sleep 1
	gprof blur_par gmon.out > outputfiles_par/im2_gprof.txt
	sleep 1
	rm gmon.out

target3:
	./blur_par 15 data/im3.ppm outputfiles_par/im3.ppm 8
	sleep 1
	gprof blur_par gmon.out > outputfiles_par/im3_gprof.txt
	sleep 1
	rm gmon.out

target4:
	./blur_par 15 data/im4.ppm outputfiles_par/im4.ppm 8
	sleep 1
	gprof blur_par gmon.out > outputfiles_par/im4_gprof.txt
	sleep 1
	rm gmon.out

blur_par: clean matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters.o -o blur_par


memtest: clean blur_valgrind memclean mem_test

blur_valgrind: matrix ppm filters blur.cpp
	$(CXX) -std=c++17 -g -Wunused -Wall -Wunused -pthread blur.cpp matrix.o ppm.o filters.o -o blur_par

mem_test:
	valgrind --tool=massif --massif-out-file=massif_im1.out ./blur_par 15 data/im1.ppm outputfiles_par/im1.ppm 8
	sleep 1
	valgrind --tool=massif --massif-out-file=massif_im2.out ./blur_par 15 data/im2.ppm outputfiles_par/im2.ppm 8
	sleep 1
	valgrind --tool=massif --massif-out-file=massif_im3.out ./blur_par 15 data/im3.ppm outputfiles_par/im3.ppm 8
	sleep 1
	valgrind --tool=massif --massif-out-file=massif_im4.out ./blur_par 15 data/im4.ppm outputfiles_par/im4.ppm 8

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
	rm -rf blur blur_par threshold *.ppm *.o *.dSYM 2> /dev/null

testclean:
	rm outputfiles_par/*

memclean:
	rm outputfiles_par/*.ppm massif*