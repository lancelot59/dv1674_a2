# Author: David Holmqvist <daae19@student.bth.se>

CXX=g++
CXXFLAGS=-std=c++17 -g -pg -Wunused -Wall -pthread

IMAGES = im1 im2 im3 im4
THREADS = 1 2 4 8 12 16 32

run: clean all seq_test par_test

run_mem_test: seq_memtest par_memtest

all: blur blur_par

threshold: matrix ppm filters threshold.cpp
	$(CXX) $(CXXFLAGS) threshold.cpp matrix.o ppm.o filters.o -o threshold

blur: matrix ppm filters blur.cpp
	$(CXX) $(CXXFLAGS) blur.cpp matrix.o ppm.o filters.o -o blur

blur_par: matrix ppm filters blur_par.cpp
	$(CXX) $(CXXFLAGS) blur_par.cpp matrix.o ppm.o filters.o -o blur_par

filters: matrix filters.hpp filters.cpp
	$(CXX) $(CXXFLAGS) -c filters.cpp -o filters.o

matrix: matrix.hpp matrix.cpp
	$(CXX) $(CXXFLAGS) -c matrix.cpp -o matrix.o

ppm: ppm.hpp ppm.cpp
	$(CXX) $(CXXFLAGS) -c ppm.cpp -o ppm.o

clean:
	rm -rf blur blur_par threshold *.ppm *.o *.dSYM 2> /dev/null

seq_test: seq_testclean
	./blur 15 data/im1.ppm outputfiles/im1.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im1_gprof.txt
	rm gmon.out
	sleep 1
	./blur 15 data/im2.ppm outputfiles/im2.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im2_gprof.txt
	rm gmon.out
	sleep 1
	./blur 15 data/im3.ppm outputfiles/im3.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im3_gprof.txt
	rm gmon.out
	sleep 1
	./blur 15 data/im4.ppm outputfiles/im4.ppm
	sleep 1
	gprof blur gmon.out > outputfiles/im4_gprof.txt
	rm gmon.out

seq_memtest: clean blur_seq_valgrind mem_test

blur_seq_valgrind: matrix ppm filters blur.cpp
	$(CXX) -std=c++17 -g -Wunused -Wall -Wunused -pthread blur.cpp matrix.o ppm.o filters.o -o blur

mem_test:
	rm outputfiles/*.out
	valgrind --tool=massif --massif-out-file=outputfiles/massif_im1.out ./blur 15 data/im1.ppm outputfiles/im1.ppm
	sleep 1
	valgrind --tool=massif --massif-out-file=outputfiles/massif_im2.out ./blur 15 data/im2.ppm outputfiles/im2.ppm
	sleep 1
	valgrind --tool=massif --massif-out-file=outputfiles/massif_im3.out ./blur 15 data/im3.ppm outputfiles/im3.ppm
	sleep 1
	valgrind --tool=massif --massif-out-file=outputfiles/massif_im4.out ./blur 15 data/im4.ppm outputfiles/im4.ppm

seq_testclean:
	rm -f outputfiles/*.txt outputfiles/*.ppm

par_test: par_testclean
	@for img in $(IMAGES); do \
		for thread in $(THREADS); do \
			echo "Running gprof test on $$img using $$thread threads"; \
			./blur_par 15 data/$$img.ppm outputfiles_par/$${img}_$${thread}.ppm $$thread; \
			sleep 1; \
			gprof blur_par gmon.out > outputfiles_par/gprof_$${img}_$${thread}.txt; \
			rm gmon.out; \
		done; \
	done

par_testclean:
	rm -f outputfiles_par/*.txt outputfiles_par/*.ppm

par_memtest: clean blur_par_valgrind par_mem_test

blur_par_valgrind: matrix ppm filters blur.cpp
	$(CXX) -std=c++17 -g -Wunused -Wall -Wunused -pthread blur_par.cpp matrix.o ppm.o filters.o -o blur_par

par_mem_test:
	rm outputfiles_par/*.out
	@for img in $(IMAGES); do \
		for thread in $(THREADS); do \
			echo "Running valgrind test on $$img using $$thread threads"; \
			valgrind \
				--tool=massif \
				--massif-out-file=outputfiles_par/massif_$${img}_$${thread}.out \
				./blur_par 15 data/$$img.ppm outputfiles_par/$${img}_$${thread}.ppm $$thread; \
			sleep 1; \
		done; \
	done
