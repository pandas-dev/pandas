CFLAGS=-O2 -std=gnu99 -I/usr/include/python2.7 -I/usr/local/lib/python2.7/dist-packages/numpy/core/include/
OBJS=test_period.o period.o
TARGET=test_period

test_period: ${OBJS}
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lpython2.7

clean:
	rm -f ${OBJS}

test_period.o: test_period.c period.o
	$(CC) ${CFLAGS} -c $<
        
period.o: period.c
	$(CC) $(CFLAGS) -c $< -lpython2.7


