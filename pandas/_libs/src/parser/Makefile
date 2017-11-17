PYTHONBASE = /Library/Frameworks/EPD64.framework/Versions/Current
NUMPY_INC = /Library/Frameworks/EPD64.framework/Versions/7.1/lib/python2.7/site-packages/numpy/core/include
PYTHON_INC = -I$(PYTHONBASE)/include/python2.7 -I$(NUMPY_INC)
PYTHON_LINK = -L$(PYTHONBASE)/lib -lpython

SOURCES = conversions.c parser.c str_to.c

check-syntax:
	gcc -g $(PYTHON_INC) -o /dev/null -S ${CHK_SOURCES}

test: $(SOURCES)
	gcc $(PYTHON_INC) -o test $(SOURCES)
	./test