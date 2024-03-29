#
#  Copyright (C) 2018 Remi Imbach
#
#  This file is part of Ccluster.
#
#  Ccluster is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License (LGPL) as published
#  by the Free Software Foundation; either version 2.1 of the License, or
#  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
#

# PROJECT=geom
CC=gcc
# CCFLAGS = -Wall -ansi -O3  -fPIC
CCFLAGS = -Wall -g -std=c99 -O3 -fno-stack-protector
# CCFLAGS =-std=c99 -pedantic -Wall -O2 -funroll-loops -g -mpopcnt

SRC_PATH = ../src
INC_PATH = ../src
SRC_DIRS= $(SRC_PATH)

FLINT_INCLUDE_PATH = /work/Ccluster/build/include
ARB_INCLUDE_PATH = /work/libs/arb

INCLUDE_PATH = -I$(FLINT_INCLUDE_PATH) -I$(ARB_INCLUDE_PATH) -I$(INC_PATH)
LIBRARY_PATH = -L/work/libs/flint2 -L/work/libs/arb -L/work/Ccluster
LIBRARY_FLAG = -lgmp -lflint -larb -lccluster -lm


SRCS=$(foreach x, $(SRC_DIRS), $(wildcard $(x)/*.c))
OBJS=$(SRCS:%.c=%.o)
HEADERS=$(foreach x, $(INC_PATH), $(wildcard $(x)/*.h))

TEST=$(foreach x, ., $(wildcard $(x)/*.c))
TEST_OBJS=$(TEST:%.c=%.o)
TEST_EXECS=$(TEST:%.c=%)

test: test_realRat test_realRat_poly test_compRat test_compRat_poly test_realApp test_compApp test_app_rat test_compApp_poly test_app_rat_poly test_cacheApp test_compDsk test_compBox test_parallel_taylor test_parallel_convol_taylor test_taylor_shift_convolution test_gen_list test_compBox_list test_connCmp test_strategies test_chronos test_counters test_metadatas test_subdBox test_ciso_algo test_graeffe test_c0_pred test_TS test_doubRealApp test_doubCompApp test_doubCompApp_poly test_doubCompApp_poly_mul test_doubCompApp_poly_mulk test_doubCompApp_poly_sqr test_doubCompApp_poly_sqrk test_doubCompApp_poly_TS test_doubCompApp_poly_G test_powerSums test_debug test_multipoint_eval

clean :
	rm -rf $(TEST_OBJS) $(TEST_EXECS)

.c.o :
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ -c $< $(INCLUDE_PATH)
	
test_parallel: test_parallel.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lpthread
	
test_realRat: test_realRat.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_realRat_poly: test_realRat_poly.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)

test_compRat: test_compRat.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)

test_compRat_poly: test_compRat_poly.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)

test_realApp: test_realApp.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_compApp: test_compApp.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_app_rat: test_app_rat.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_compApp_poly: test_compApp_poly.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_app_rat_poly: test_app_rat_poly.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_cacheApp: test_cacheApp.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lpthread
	
test_compDsk: test_compDsk.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_compBox: test_compBox.o 
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_parallel_taylor: test_parallel_taylor.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_parallel_convol_taylor: test_parallel_convol_taylor.o
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_taylor_shift_convolution: test_taylor_shift_convolution.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_gen_list: test_gen_list.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_compBox_list: test_compBox_list.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_connCmp: test_connCmp.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_strategies: test_strategies.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_chronos: test_chronos.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_counters: test_counters.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_metadatas: test_metadatas.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_subdBox: test_subdBox.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_ciso_algo: test_ciso_algo.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_graeffe: test_graeffe.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_c0_pred: test_c0_pred.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_TS: test_TS.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG)
	
test_doubRealApp: test_doubRealApp.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp: test_doubCompApp.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly: test_doubCompApp_poly.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_mul: test_doubCompApp_poly_mul.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_mulk: test_doubCompApp_poly_mulk.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_sqr: test_doubCompApp_poly_sqr.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_sqrk: test_doubCompApp_poly_sqrk.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_TS: test_doubCompApp_poly_TS.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_doubCompApp_poly_G: test_doubCompApp_poly_G.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
testSpiralDac: testSpiralDac.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_powerSums: test_powerSums.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm

test_debug: test_debug.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm

test_multipoint_eval: test_multipoint_eval.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_fpri: test_fpri.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_fpri_run: test_fpri
	@echo "Testing fpri"
	test_fpri
	valgrind --tool=memcheck --leak-check=yes test_fpri
	
test_fpci: test_fpci.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
test_fpci_run: test_fpci
	@echo "Testing fpci"
	test_fpci
# 	valgrind --tool=memcheck --leak-check=yes test_fpci
	
test_fpri_poly: test_fpri_poly.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
deg=100
test_fpri_poly_run: test_fpri_poly
	@echo "Testing fpri_poly"
# 	echo "7 49/128 -7/16 -213/128 35/16 9/32 -7/4 1" > test_poly.txt
	../bin/genPolFile randomDense ${deg} test_poly.txt
	test_fpri_poly test_poly.txt
# 	valgrind --tool=memcheck --leak-check=yes test_fpri_poly test_poly.txt
	rm -rf test_poly.txt
	
test_fpci_poly: test_fpci_poly.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
deg=100
nbt=10
test_fpci_poly_run: test_fpci_poly
	@echo "Testing fpci_poly"
# 	echo "7 49/128 -7/16 -213/128 35/16 9/32 -7/4 1" > test_poly.txt
	../bin/genPolFile randomDense ${deg} test_poly.txt
# 	../bin/genPolFile randomSparse ${deg} test_poly.txt -n ${nbt}
	test_fpci_poly test_poly.txt
# 	valgrind --tool=memcheck --leak-check=yes test_fpri_poly test_poly.txt
	rm -rf test_poly.txt

time_fpci_sparse_eval: time_fpci_sparse_eval.c
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $< $(INCLUDE_PATH) $(LIBRARY_PATH) $(LIBRARY_FLAG) -lm
	
time_fpci_sparse_eval_run: time_fpci_sparse_eval
	@echo "Testing fpci_poly"
# 	echo "7 49/128 -7/16 -213/128 35/16 9/32 -7/4 1" > test_poly.txt
	../bin/genPolFile randomDense ${deg} test_poly.txt
# 	../bin/genPolFile randomSparse ${deg} test_poly.txt -n ${nbt}
	time_fpci_sparse_eval test_poly.txt
# 	valgrind --tool=memcheck --leak-check=yes test_fpri_poly test_poly.txt
	rm -rf test_poly.txt
