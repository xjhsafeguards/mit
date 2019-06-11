
TOPDIR=$(shell git rev-parse --show-toplevel)/
FOLDERS=Cell Math Physics Utilit #old-Cell 

# Compiler
CC=g++

# Define and include flags
DFLAGS= -DNDEBUG #-DQUANT 
IFLAGS=-I$(TOPDIR)/eigen-eigen-323c052e1731/ $(patsubst %,-I$(TOPDIR)/mit/%,$(FOLDERS)) #-I$(TOPDIR)/mit/include/ 

# C flags
CFLAGS=-O3 -std=c++11 $(DFLAGS) $(IFLAGS) 

# Linker
LD=g++
LDFLAGS=#-static
LD_LIBS=

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv
#AR             = ar
#ARFLAGS        = ruv
