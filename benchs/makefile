#
#  Copyright (C) 2018 Remi Imbach, Chee K. Yap
#
#  This file is part of Ccluster.
#
#  Ccluster is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License (LGPL) as published
#  by the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
#

#Makefile

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# VARIABLES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# TESTING VARIABLES

exe=
exe=.exe

# for benchMignotte:
degree=64
bitsize=14
box=0,1,0,1,20,1 		# (x,y), width
epsilon=1,200			# 1/100
strategy=15
verbosity=2



# COMPILING VARIABLES
CC=gcc
CCFLAGS = -Wall -std=c99 -O2
# CCFLAGS = -Wall -g

SRC_PATH = ../src
INC_PATH = ../src
SRC_DIRS= $(SRC_PATH)

FLINT_INCLUDE_PATH = /usr/local/include/flint
ARB_INCLUDE_PATH = /usr/local/include

INCLUDE_PATH = -I$(FLINT_INCLUDE_PATH) -I$(ARB_INCLUDE_PATH) -I$(INC_PATH)
#LIBRARY_PATH = -L /usr/local/lib -L ..
LIBRARY_PATH = -L /usr/local/lib:..:.
LIBRARY_FLAG = -lflint -larb -lccluster -lm


SRCS=$(foreach x, $(SRC_DIRS), $(wildcard $(x)/*.c))
OBJS=$(SRCS:%.c=%.o)
HEADERS=$(foreach x, $(INC_PATH), $(wildcard $(x)/*.h))

TEST=$(foreach x, ., $(wildcard $(x)/*.c))
TEST_OBJS=$(TEST:%.c=%.o)
TEST_EXECS=$(TEST:%.c=%)

BENCH_DIRS= ./Bernoulli\
            ./Cluster\
            ./Mignotte\
            ./MignotteMul\
            ./Spiral\
            ./Wilkinson\
            ./WilkMul\
            ./WilkRat
BENCH_DIRS_FILES=$(foreach x, $(BENCH_DIRS), $(wildcard $(x)/*))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TARGETS:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default: test

#==================================================
# TESTING 
#==================================================
testMignotte1 tM1:
	./benchMignotte.exe 64 14 0,1,0,1,20,1 1,2000 15 2 

testMignotte tM:
	./benchMignotte${exe} \
		${degree} \
		${bitsize} \
		${box} \
		${epsilon} \
		${strategy} \
		${verbosity} \

testMignotte2 tM2:
	make testMignotte \
	    bitsize=200 \
	    epsilon=1,1000000000000 \
	    box=0,1,0,1,2,1 \
	    verbosity=3

testBernoulli tB:
	./benchBernoulli${exe} \
		${degree} \
		${box} \
		${epsilon} \
		${strategy} \
		${verbosity} \

testBernoulli1 tB1:
	./benchBernoulli.exe 64 0,1,0,1,20,1  1,2000 15 2 

testWilkinson tW:
	./benchWilkinson${exe} \
		${degree} \
		${box} \
		${epsilon} \
		${strategy} \
		${verbosity} \

testWilkinson1 tW1:
	./benchWilkinson.exe 64 0,1,0,1,20,1  1,2000 15 2 

#==================================================
# TABLES 
#==================================================

tableBernoullim:
	./tableBernoulli.sh --purge --degrees="64 128 191 256 383 512 767" --stop-when-compact --anticipate 2>/dev/null > tableBernoulli.txt &
	
tableV1V2V3m:
	./tableV1V2V3transposed.sh --purge --degree=64 --bitsize=14 --nbSols=11 --pow=3 --itts=4 2>/dev/null
	
tableV3V4m:
	./tableV3V4transposed.sh --purge --degree="64 128 256" --bitsize=14 --nbSols="11 12 13" --pow=3 --itts="3 4 5" 2>/dev/null
	
tableV4V4escm:
	./tableV4V4esc.sh --degree=64 --bitsize=14 --epsilon=-53 --eps1=-53 --eps2=-530 --eps3=-5300  --nbSols=11 --pow=3 --itts=4 2>/dev/null

#==================================================
# COMPILING 
#==================================================
benchs: benchMignotte benchMignotteMul \
    benchBernoulli benchWilkinson \
    benchWilkRat benchSpiral \
    benchCluster benchWilkMul \
    benchMandelbrot benchLaguerre\
    gen

gen: genPolFile

#==================================================
# HOUSEKEEPING 
#==================================================
# now in ../Makefile with make benchs
# print: 
# 	@echo $(BENCH_DIRS_FILES)
# 
# clean :
# 	rm -rf $(TEST_OBJS) $(TEST_EXECS)
# 	
# cleanBenchs : 
# 	rm -rf $(BENCH_DIRS_FILES)
# 
# .c.o :
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ -c $< $(INCLUDE_PATH)
# 	
# benchMignotte: benchMignotte.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) \
# 	    	-o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG) \
# 	
# debug: benchMignotte.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG) \
# 	    	-o $@ $< $(INCLUDE_PATH) \
# 	
# benchMignotteMul: benchMignotteMul.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchWilkinson: benchWilkinson.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchBernoulli: benchBernoulli.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchWilkRat: benchWilkRat.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchSpiral: benchSpiral.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchCluster: benchCluster.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# benchWilkMul: benchWilkMul.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
# 	
# genPolFile: genPolFile.o
# 	@echo "Building $@"
# 	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
# 	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)

#==================================================
# END 
#==================================================
