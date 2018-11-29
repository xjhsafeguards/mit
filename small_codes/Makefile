CC=g++
OPTION=-O3 -std=c++11 

SOURCES=${wildcard */*.cpp}
OBJS=${SOURCES:.cpp=.o}

all:${OBJS}

.cpp.o:
	${CC} ${OPTION} -c $< -o $@

clean:
	@-rm -f *.o *.x *.log
	@echo "Object cleaned!"
	@cd old && make clean
