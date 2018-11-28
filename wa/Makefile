CC=g++
OPTION=-O3 -std=c++11 -DN_TEST 

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=wa.x

all:${EXECUTABLE}

.cpp.o:
	${CC} ${OPTION} -c $< -o $@

${EXECUTABLE}:${OBJS}
	${CC} ${OPTION} -o ../../${EXECUTABLE} ${OBJS}

.PHONY : old
old:
	@cd old && $(MAKE)

.PHONY : clean
clean:
	@-rm -f *.o *.x *.log
	@echo "Object cleaned!"
	@cd old && make clean
