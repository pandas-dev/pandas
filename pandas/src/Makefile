CFLAGS=-O2 -g -std=gnu99 $(shell python-config --include) -I$(shell python -c "import numpy; print numpy.get_include()")
LDFLAGS=-lpython2.7
OBJS=test_period.o period.o
TARGET=test_period

test_period: ${OBJS}
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS) 

clean:
	rm -f ${OBJS}

test_period.o: test_period.c period.o
	$(CC) ${CFLAGS} -c $<
        
period.o: period.c
	$(CC) $(CFLAGS) -c $<


