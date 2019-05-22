include make.inc

FOLDERS=Utilit Math

SOURCES=${wildcard *.cpp}
OBJS=${SOURCES:.cpp=.o}
EXECUTABLE=mit.x

all:${EXECUTABLE}

.cpp.o:
	${CC} ${CFLAGS} -c $< -o $@ 

${EXECUTABLE}:${OBJS} lib
	$(LD) $(LDFLAGS) $< ${wildcard lib/*.a} -o $@ 

%.x: %.o lib
	$(LD) $(LDFLAGS) $< ${wildcard lib/*.a} -o $@ 

.PHONY : lib
lib:
	@for dir in ${FOLDERS}\
	; do \
		if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= all ) \
		fi \
	done

.PHONY : clean
clean:
	@for dir in ${FOLDERS}\
	; do \
		if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
		fi \
	done
	@-rm -f *.x *.o *.log lib/*.a
	@echo "All object cleaned!"
