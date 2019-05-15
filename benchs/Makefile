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

# COMPILING VARIABLES
CC=gcc
CCFLAGS = -Wall -std=c99 -O2
# CCFLAGS = -Wall -g
QUIET_CC  = @echo '   ' CC  ' ' $@;

SRC_PATH = ../src
INC_PATH = ../src
SRC_DIRS= $(SRC_PATH)

#FLINT_INCLUDE_PATH = /usr/local/include/flint
#ARB_INCLUDE_PATH = /usr/local/include
FLINT_INCLUDE_PATH = /home/imbach/Ccluster/build/include/
ARB_INCLUDE_PATH = /home/imbach/libs/arb

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

BENCH_DIRS= ./tableBernoulli\
	    ./tableLaguerre\
	    ./tableMandelbrot\
	    ./tableMignotte\
	    ./tableMignotteGen\
	    ./tableNestedClusters\
	    ./tableSpiral\
	    ./tableV1V2V3\
	    ./tableV3V4\
	    ./tableV4esc_epsilon\
	    ./tableV4V4esc\
	    ./tableWilkinson\
	    ./tableWilkMul\
	    ./tableWilkRat
BENCH_DIRS_FILES=$(foreach x, $(BENCH_DIRS), $(wildcard $(x)/*))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TARGETS:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default: build

#==================================================
# TABLES 
#==================================================

tableBernoulliM:
	./tableBernoulli.sh --purge --degrees="64 128 191 256 383 512 767"\
	                    --stop-when-compact --anticipate --mpsolve=0\
	                    2>/dev/null
	
tableV1V2V3M:
	./tableV1V2V3.sh --purge --degree=64\
	--bitsize=14 --nbSols=11 --pow=3 --itts=4\ 
	2>/dev/null
	
tableV3V4M:
	./tableV3V4.sh --purge --degree="64 128 256 512" \
	--bitsize=14 --nbSols="11 12 13 14" --pow=3 --itts="3 4 5"\
	2>/dev/null
	
tableV4V4escM:
	./tableV4V4esc.sh --degree=64\
	--bitsize=14 --epsilon=-53 --eps1=-53 --eps2=-530 --eps3=-5300\
	--nbSols=11 --pow=3 --itts=4\
	2>/dev/null

tableICMS: tableBernoulliM\
           tableV1V2V3M\
           tableV3V4M\
           tableV4V4escM
           
           
#==================================================
# COMPILING 
#==================================================
build: genPolFile

#==================================================
# HOUSEKEEPING 
#==================================================
print: 
	@echo $(BENCH_DIRS_FILES)
# 
clean :
	rm -rf $(TEST_OBJS) $(TEST_EXECS)
# 	
cleanBenchs : 
	rm -rf $(BENCH_DIRS)

.c.o :
# 	@echo "Building $@"
	$(QUIET_CC) $(CC) $(CCFLAGS) -o $@ -c $< $(INCLUDE_PATH)
	
	
genPolFile: genPolFile.o
# 	@echo "Building $@"
	$(QUIET_CC) $(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) \
	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)

#==================================================
# END 
#==================================================
