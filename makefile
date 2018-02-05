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

PROJECT=ccluster
CC=gcc
# CCFLAGS = -Wall -ansi -O3  -fPIC
CCFLAGS = -Wall -fPIC -std=c99 -O2
# Chee:
#CCFLAGS = -Wall -std=c99 -O2
# CCFLAGS = -Wall -fPIC -g

SRC_PATH = ./src
INC_PATH = ./src
          
SRC_DIRS= $(SRC_PATH)\
          $(SRC_PATH)/base\
          $(SRC_PATH)/numbers\
          $(SRC_PATH)/polynomials\
          $(SRC_PATH)/tstar\
          $(SRC_PATH)/caches\
          $(SRC_PATH)/geometry\
          $(SRC_PATH)/lists\
          $(SRC_PATH)/metadatas\
          $(SRC_PATH)/newton\
          $(SRC_PATH)/ccluster
          
BUILD_DIR = ./build/

FLINT_INCLUDE_PATH = /usr/local/include/flint
ARB_INCLUDE_PATH = /usr/local/include/

INCLUDE_PATH = -I$(FLINT_INCLUDE_PATH) -I$(ARB_INCLUDE_PATH) -I$(INC_PATH)
LIBRARY_PATH = -L /usr/local/lib
LIBRARY_FLAG = -lflint -larb -lpthread

SRCS=$(foreach x, $(SRC_DIRS), $(wildcard $(x)/*.c))
OBJS=$(SRCS:%.c=%.o)
HEADERS=$(foreach x, $(INC_PATH), $(wildcard $(x)/*.h))

print: 
	@echo $(OBJS)
	@echo $(SRCS)

DYNAMIC= -shared
DYNAMIC_LIB_LINUX= lib$(PROJECT).so
DYNAMIC_LIB_WINDOWS= lib$(PROJECT).dll

all : $(DYNAMIC_LIB_WINDOWS)
allinux: $(DYNAMIC_LIB_LINUX)

$(DYNAMIC_LIB_WINDOWS) : $(OBJS)
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $(DYNAMIC) $(OBJS) \
	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)
	    	
$(DYNAMIC_LIB_LINUX) : $(OBJS)
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ $(DYNAMIC) $(OBJS) \
	    	$(LIBRARY_PATH) $(LIBRARY_FLAG)

.c.o :
	@echo "Building $@"
	$(CC) $(CCFLAGS) -o $@ -c $< $(INCLUDE_PATH)
# 	mv $@ $(BUILD_DIR)
	
clean: 
	@echo "Cleaning all"
	rm -rf $(OBJS) $(DYNAMIC_LIB_LINUX) $(DYNAMIC_LIB_WINDOWS)
	
	
