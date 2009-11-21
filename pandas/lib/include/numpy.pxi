# :Author:    Travis Oliphant

cdef extern from "numpy/arrayobject.h":

	cdef enum NPY_TYPES:
		NPY_BOOL
		NPY_BYTE
		NPY_UBYTE
		NPY_SHORT
		NPY_USHORT
		NPY_INT
		NPY_UINT
		NPY_LONG
		NPY_ULONG
		NPY_LONGLONG
		NPY_ULONGLONG
		NPY_FLOAT
		NPY_DOUBLE
		NPY_LONGDOUBLE
		NPY_CFLOAT
		NPY_CDOUBLE
		NPY_CLONGDOUBLE
		NPY_OBJECT
		NPY_STRING
		NPY_UNICODE
		NPY_VOID
		NPY_NTYPES
		NPY_NOTYPE

	cdef enum requirements:
		NPY_CONTIGUOUS
		NPY_FORTRAN
		NPY_OWNDATA
		NPY_FORCECAST
		NPY_ENSURECOPY
		NPY_ENSUREARRAY
		NPY_ELEMENTSTRIDES
		NPY_ALIGNED
		NPY_NOTSWAPPED
		NPY_WRITEABLE
		NPY_UPDATEIFCOPY
		NPY_ARR_HAS_DESCR

		NPY_BEHAVED
		NPY_BEHAVED_NS
		NPY_CARRAY
		NPY_CARRAY_RO
		NPY_FARRAY
		NPY_FARRAY_RO
		NPY_DEFAULT

		NPY_IN_ARRAY
		NPY_OUT_ARRAY
		NPY_INOUT_ARRAY
		NPY_IN_FARRAY
		NPY_OUT_FARRAY
		NPY_INOUT_FARRAY

		NPY_UPDATE_ALL

	cdef enum NPY_SORTKIND:
		NPY_QUICKSORT
		NPY_HEAPSORT
		NPY_MERGESORT

	cdef enum defines:
		NPY_MAXDIMS

	ctypedef struct npy_cdouble:
		double real
		double imag

	ctypedef struct npy_cfloat:
		double real
		double imag

	ctypedef int npy_intp

	ctypedef extern class numpy.dtype [object PyArray_Descr]:
		cdef int type_num, elsize, alignment
		cdef char type, kind, byteorder, hasobject
		cdef object fields, typeobj
		cdef object f

	ctypedef extern class numpy.ndarray [object PyArrayObject]:
		cdef char *data
		cdef int nd
		cdef npy_intp *dimensions
		cdef npy_intp *strides
		cdef object base
		cdef dtype descr
		cdef int flags

	ctypedef extern class numpy.flatiter [object PyArrayIterObject]:
		cdef int  nd_m1
		cdef npy_intp index, size
		cdef ndarray ao
		cdef char *dataptr

	ctypedef extern class numpy.broadcast [object PyArrayMultiIterObject]:
		cdef int numiter
		cdef npy_intp size, index
		cdef int nd
		cdef npy_intp *dimensions
		cdef void **iters

	object PyArray_ZEROS(int ndims, npy_intp* dims, int type_num, int fortran)
	object PyArray_EMPTY(int ndims, npy_intp* dims, int type_num, int fortran)
	dtype PyArray_DescrFromTypeNum(NPY_TYPES type_num)
	object PyArray_SimpleNew(int ndims, npy_intp* dims, NPY_TYPES type_num)
	int PyArray_Check(object obj)
	object PyArray_ContiguousFromAny(object obj, int type,
		int mindim, int maxdim)
	object PyArray_ContiguousFromObject(object obj, NPY_TYPES type,
		int mindim, int maxdim)
	npy_intp PyArray_SIZE(ndarray arr)
	npy_intp PyArray_NBYTES(ndarray arr)
	npy_intp *PyArray_DIMS(object arr)
	void *PyArray_DATA(ndarray arr)
	int PyArray_ITEMSIZE(object arr)

	int PyArray_ISBEHAVED(ndarray arr)

	int PyArray_TYPE(object arr)

	object PyArray_GETITEM(object obj, void* itemptr)
	int PyArray_SETITEM(object obj, void* itemptr, object item)
	object PyArray_FromAny(object obj, dtype newtype, int mindim, int maxdim,
						   int requirements, object context)
	object PyArray_FROMANY(object obj, NPY_TYPES type_num, int min,
						   int max, int requirements)
	object PyArray_NewFromDescr(object subtype, dtype newtype, int nd,
								npy_intp* dims, npy_intp* strides, void* data,
								int flags, object parent)

	object PyArray_FROM_OTF(object obj, NPY_TYPES type, int flags)
	object PyArray_EnsureArray(object)

	object PyArray_MultiIterNew(int n, ...)

	char *PyArray_MultiIter_DATA(broadcast multi, int i)
	void PyArray_MultiIter_NEXTi(broadcast multi, int i)
	void PyArray_MultiIter_NEXT(broadcast multi)

	object PyArray_IterNew(object arr)
	void PyArray_ITER_NEXT(flatiter it)
	void PyArray_ITER_GOTO(flatiter it, npy_intp* i)
	void PyArray_ITER_GOTO1D(flatiter it, npy_intp i)

	void import_array()
