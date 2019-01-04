CC=g++
OPTION=-O3 -std=c++11 

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=mit.x

all:${EXECUTABLE}

.cpp.o:
	${CC} ${OPTION} -c $< -o $@

${EXECUTABLE}:${OBJS}
	${CC} ${OPTION} -o ../../${EXECUTABLE} ${OBJS}

.PHONY : clean
clean:
	@-rm -f *.o *.x *.log
	@echo "Object cleaned!"
