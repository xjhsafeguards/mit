CC=g++
OPTION=-O3 -std=c++11 -I ../eigen-eigen-323c052e1731/ 

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=mit.x

all:${EXECUTABLE}

.cpp.o:
	${CC} ${OPTION} -c $< -o $@

${EXECUTABLE}:${OBJS}
	${CC} ${OPTION} -o ${EXECUTABLE} ${OBJS}

.PHONY : clean
clean:
	@-rm -f *.o *.x *.log
	@echo "Object cleaned!"
