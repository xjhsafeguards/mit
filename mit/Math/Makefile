#Makefile for Math

include ../make.inc

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=../lib/mymath.a

all:${EXECUTABLE}

.cpp.o:
	${CC} ${CFLAGS} -c $< -o $@

${EXECUTABLE}:${OBJS}
	${AR} ${ARFLAGS} ${EXECUTABLE} ${OBJS}

.PHONY : clean

clean:
	@-rm -f *.o *.a *.log
	@echo "Math object cleaned!"
