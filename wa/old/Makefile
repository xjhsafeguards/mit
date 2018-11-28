CC=g++
OPTION=-O3 -std=c++11 

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=wao.x

all:${EXECUTABLE}

.cpp.o:
	${CC} ${OPTION} -c $< -o $@

${EXECUTABLE}:${OBJS}
	${CC} ${OPTION} -o ../../../${EXECUTABLE} ${OBJS}

clean:
	@-rm -f *.o *.x *.log
