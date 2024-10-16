# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -pg -Wunused -Wall -Wunused -mavx -pthread

all: pearson

run: clean all test

test: target1 target2 target3 target4

target1:
	./pearson data/128.data outputfiles/128.data
	sleep 1
	gprof ./pearson gmon.out > outputfiles/128gprof.txt
	sleep 1
	rm gmon.out

target2: target1
	./pearson data/256.data outputfiles/256.data
	sleep 1 
	gprof ./pearson gmon.out > outputfiles/256gprof.txt
	sleep 1
	rm gmon.out

target3: target2
	./pearson data/512.data outputfiles/512.data
	sleep 1
	gprof ./pearson gmon.out > outputfiles/512gprof.txt
	sleep 1
	rm gmon.out

target4:target3
	./pearson data/1024.data outputfiles/1024.data
	sleep 1
	gprof ./pearson gmon.out > outputfiles/1024gprof.txt
	sleep 1
	rm gmon.out
	
debugg:
	gdb ./pearson
	data/128.data outputfiles/debugg128.data

memtest: clean all
	rm outputfiles/*.data massif*
	valgrind --tool=massif ./pearson data/128.data outputfiles/128.data
	valgrind --tool=massif ./pearson data/256.data outputfiles/256.data
	valgrind --tool=massif ./pearson data/512.data outputfiles/512.data
	valgrind --tool=massif ./pearson data/1024.data outputfiles/1024.data



pearson: vector dataset analysis pearson.cpp
	$(CXX) $(CXXFLAGS) pearson.cpp vector.o dataset.o analysis.o -o pearson

analysis: vector analysis.hpp analysis.cpp
	$(CXX) $(CXXFLAGS) -c analysis.cpp -o analysis.o

dataset: vector dataset.hpp dataset.cpp
	$(CXX) $(CXXFLAGS) -c dataset.cpp -o dataset.o

vector: vector.hpp vector.cpp
	$(CXX) $(CXXFLAGS) -c vector.cpp -o vector.o

clean:
	rm -rf pearson *.o *.dSYM 2> /dev/null


testclean:
	rm outputfiles/*