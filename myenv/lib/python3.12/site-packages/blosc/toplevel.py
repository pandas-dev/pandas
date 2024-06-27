########################################################################
#
#       Created: September 22, 2010
#       Author:  The Blosc development team - blosc@blosc.org
#
########################################################################

import os
import pickle
import subprocess
import sys
from ._version import LooseVersion

from blosc import blosc_extension as _ext
import blosc


def detect_number_of_cores():
    """
    detect_number_of_cores()

    Detect the number of cores in this system.

    Returns
    -------
    out : int
        The number of cores in this system.

    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else:  # OSX:
            completed = subprocess.run(("sysctl", "-n", "hw.ncpu"),
                                       stdout=subprocess.PIPE)
            if completed.returncode == 0:
                ncpus = int(completed.stdout)
                if ncpus > 0:
                    return ncpus
    # Windows:
    elif "NUMBER_OF_PROCESSORS" in os.environ:
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    return 1  # Default


def set_nthreads(nthreads):
    """
    set_nthreads(nthreads)

    Set the number of threads to be used during Blosc operation.

    Parameters
    ----------
    nthreads : int
        The number of threads to be used during Blosc operation.

    Returns
    -------
    out : int
        The previous number of used threads.

    Raises
    ------
    ValueError
        If nthreads is larger that the maximum number of threads blosc can use.

    Notes
    -----
    The number of threads for Blosc is the maximum number of cores
    detected on your machine (via `detect_number_of_cores`).  In some
    cases Blosc gets better results if you set the number of threads
    to a value slightly below than your number of cores.

    Examples
    --------
    Set the number of threads to 2 and then to 1:

    >>> oldn = blosc.set_nthreads(2)
    >>> blosc.set_nthreads(1)
    2

    """
    if nthreads > blosc.MAX_THREADS:
        raise ValueError("the number of threads cannot be larger than %d" %
                         blosc.MAX_THREADS)

    blosc.nthreads = nthreads
    return _ext.set_nthreads(nthreads)


def set_blocksize(blocksize):
    """set_blocksize(blocksize)

    Force the use of a specific blocksize.  If 0, an automatic
    blocksize will be used (the default).

    Notes
    -----

    This is a low-level function and is recommended for expert users only.

    Examples
    --------

    >>> blosc.set_blocksize(512)
    >>> blosc.set_blocksize(0)

    """

    _ext.set_blocksize(blocksize)


def get_blocksize():
    """
    get_blocksize()

    Get the blocksize currently used.

    Examples
    --------

    >>> blosc.set_blocksize(0)
    >>> blosc.get_blocksize()
    0

    """
    return _ext.get_blocksize()


def set_releasegil(gilstate):
    """
    set_releasegil( gitstate )

    Sets a boolean on whether to release the Python global inter-lock (GIL)
    during c-blosc compress and decompress operations or not.  This defaults
    to False.

    Notes
    -----

    Designed to be used with larger chunk sizes and a ThreadPool.  There is a
    small performance penalty with releasing the GIL that will more harshly
    penalize small block sizes.

    Examples
    --------

    >>> oldReleaseState = blosc.set_releasegil(True)

    """
    gilstate = bool(gilstate)
    return _ext.set_releasegil(gilstate)


def compressor_list():
    """
    compressor_list()

    Returns a list of compressors available in C library.

    Parameters
    ----------
    None

    Returns
    -------
    out : list
        The list of names.
    """
    return _ext.compressor_list().split(',')


def code_to_name(code):
    """
    code_to_name(code)

    Return the compressor name of a compressor code.

    Parameters
    ----------
    code : int
        The compressor code.

    Returns
    -------
    out : str
        The compressor name.
    """
    return _ext.code_to_name(code)


def name_to_code(name):
    """
    name_to_code(name)

    Return the compressor code of a compressor name.

    Parameters
    ----------
    name : str
        The compressor name.

    Returns
    -------
    out : int
        The compressor code.
    """
    return _ext.name_to_code(name)


def clib_info(cname):
    """
    clib_info(cname)

    Return info for compression libraries in C library.

    Parameters
    ----------
    cname : str
        The compressor name.

    Returns
    -------
    out : tuple
        The associated library name and version.
    """
    return _ext.clib_info(cname)


def get_clib(bytesobj):
    """
    get_clib(bytesobj)

    Return the name of the compression library for Blosc `bytesobj` buffer.

    Parameters
    ----------
    bytesobj : str / bytes
        The compressed buffer.

    Returns
    -------
    out : str
        The name of the compression library.
    """
    _check_bytesobj(bytesobj)

    return _ext.get_clib(bytesobj)


def get_cbuffer_sizes(bytesobj):
    """
    get_cbuffer_sizes(bytesobj)

    Return information about a compressed buffer:
    (number of uncompressed bytes, number of compressed bytes, blocksize).

    Parameters
    ----------
    bytesobj : str / bytes
        The compressed buffer.

    Returns
    -------
    out : tuple
        The associated uncompressed bytes, compressed bytes and blocksize

    Examples
    --------

    >>> s = b'0123456789'
    >>> c = blosc.compress(s, typesize=1)
    >>> blosc.get_cbuffer_sizes(c)
    (10, 26, 10)
    """

    return _ext.get_cbuffer_sizes(bytesobj)


def cbuffer_validate(bytesobj):
    """
    validate_cbuffer(bytesobj)

    Validate the cbuffer. Check that the cbuffer is safe to compress. Note
    that this does not guarantee that the blosc will be able to decompress the
    buffer successfully, only that it is safe to attempt to do so.

    Parameters
    ----------
    bytesobj : str / bytes
        The compressed buffer.

    Returns
    -------
    result : boolean
        True if safe, False if not.

    Examples
    --------

    >>> s = b'0123456789'
    >>> c = blosc.compress(s, typesize=1)
    >>> blosc.cbuffer_validate(c)
    True
    """

    return _ext.cbuffer_validate(bytesobj)


def free_resources():
    """
    free_resources()

    Free possible memory temporaries and thread resources.

    Returns
    -------
        out : None

    Notes
    -----
    Blosc maintain a pool of threads waiting for work as well as some
    temporary space.  You can use this function to release these
    resources when you are not going to use Blosc for a long while.

    Examples
    --------

    >>> blosc.free_resources()
    >>>
    """
    _ext.free_resources()


def _check_shuffle(shuffle):
    if shuffle not in (blosc.NOSHUFFLE, blosc.SHUFFLE, blosc.BITSHUFFLE):
        raise ValueError("shuffle can only be one of NOSHUFFLE, SHUFFLE"
                         " and BITSHUFFLE.")
    if (shuffle == blosc.BITSHUFFLE and
        LooseVersion(blosc.blosclib_version) < LooseVersion("1.8.0")):
        raise ValueError("You need C-Blosc 1.8.0 or higher for using"
                         " BITSHUFFLE.")


def _check_clevel(clevel):
    if not 0 <= clevel <= 9:
        raise ValueError("clevel can only be in the 0-9 range.")


def _check_cname(cname):
    list_cnames = compressor_list()
    if cname not in list_cnames:
        raise ValueError("cname can only be one of: %s, not '%s'" %
                         (list_cnames, cname))


def _check_typesize(typesize):
    if not 1 <= typesize <= blosc.MAX_TYPESIZE:
        raise ValueError("typesize can only be in the 1-%d range." %
                         blosc.MAX_TYPESIZE)


def _check_bytesobj(bytesobj):
    if not isinstance(bytesobj, bytes):
        raise TypeError("only string (2.x) or bytes (3.x) objects "
                        "supported as input")


def _check_byteslike(bytes_like):
    try:
        memoryview(bytes_like)
    except Exception:
        raise TypeError("Input type %s must be a bytes-like object that supports Python Buffer Protocol" % type(bytes_like))


def _check_input_length(input_name, input_len):
    if input_len > blosc.MAX_BUFFERSIZE:
        raise ValueError("%s cannot be larger than %d bytes" %
                         (input_name, blosc.MAX_BUFFERSIZE))


def _check_address(address):
    if not isinstance(address, int):
        raise TypeError("only int or long objects are supported as address")


def compress(bytesobj, typesize=8, clevel=9, shuffle=blosc.SHUFFLE,
             cname='blosclz'):
    """compress(bytesobj[, typesize=8, clevel=9, shuffle=blosc.SHUFFLE, cname='blosclz']])

    Compress bytesobj, with a given type size.

    Parameters
    ----------
    bytesobj : bytes-like object (supporting the buffer interface)
        The data to be compressed.
    typesize : int
        The data type size.
    clevel : int (optional)
        The compression level from 0 (no compression) to 9
        (maximum compression).  The default is 9.
    shuffle : int (optional)
        The shuffle filter to be activated.  Allowed values are
        blosc.NOSHUFFLE, blosc.SHUFFLE and blosc.BITSHUFFLE.  The
        default is blosc.SHUFFLE.
    cname : string (optional)
        The name of the compressor used internally in Blosc. It can be
        any of the supported by Blosc ('blosclz', 'lz4', 'lz4hc',
        'snappy', 'zlib', 'zstd' and maybe others too). The default is
        'blosclz'.

    Returns
    -------
    out : str / bytes
        The compressed data in form of a Python str / bytes object.

    Raises
    ------
    TypeError
        If bytesobj doesn't support the buffer interface.
    ValueError
        If bytesobj is too long.
        If typesize is not within the allowed range.
        If clevel is not within the allowed range.
        If cname is not a valid codec.

    Examples
    --------

    >>> import array, sys
    >>> a = array.array('i', range(1000*1000))
    >>> a_bytesobj = a.tobytes()
    >>> c_bytesobj = blosc.compress(a_bytesobj, typesize=4)
    >>> len(c_bytesobj) < len(a_bytesobj)
    True

    """

    _check_input_length('bytesobj', len(bytesobj))
    _check_typesize(typesize)
    _check_shuffle(shuffle)
    _check_clevel(clevel)
    _check_cname(cname)

    return _ext.compress(bytesobj, typesize, clevel, shuffle, cname)


def compress_ptr(address, items, typesize=8, clevel=9, shuffle=blosc.SHUFFLE,
                 cname='blosclz'):
    """compress_ptr(address, items[, typesize=8, clevel=9, shuffle=blosc.SHUFFLE, cname='blosclz']])

    Compress the data at address with given items and typesize.

    Parameters
    ----------
    address : int or long
        the pointer to the data to be compressed
    items : int
        The number of items (of typesize) to be compressed.
    typesize : int
        The data type size.
    clevel : int (optional)
        The compression level from 0 (no compression) to 9
        (maximum compression).  The default is 9.
    shuffle : int (optional)
        The shuffle filter to be activated.  Allowed values are
        blosc.NOSHUFFLE, blosc.SHUFFLE and blosc.BITSHUFFLE.  The
        default is blosc.SHUFFLE.
    cname : string (optional)
        The name of the compressor used internally in Blosc. It can be
        any of the supported by Blosc ('blosclz', 'lz4', 'lz4hc',
        'snappy', 'zlib', 'zstd' and maybe others too). The default is
        'blosclz'.

    Returns
    -------
    out : str / bytes
        The compressed data in form of a Python str / bytes object.

    Raises
    ------
    TypeError
        If address is not of type int or long.
    ValueError
        If items * typesize is larger than the maximum allowed buffer size.
        If typesize is not within the allowed range.
        If clevel is not within the allowed range.
        If cname is not within the supported compressors.

    Notes
    -----
    This function can be used anywhere that a memory address is available in
    Python. For example the Numpy "__array_interface__['data'][0]" construct,
    or when using the ctypes modules.

    Importantly, the user is responsible for making sure that the memory
    address is valid and that the memory pointed to is contiguous. Passing a
    non-valid address has a high likelihood of crashing the interpreter by
    segfault.

    Examples
    --------

    >>> import numpy
    >>> items = 7
    >>> np_array = numpy.arange(items)
    >>> c = blosc.compress_ptr(np_array.__array_interface__['data'][0], \
        items, np_array.dtype.itemsize)
    >>> d = blosc.decompress(c)
    >>> np_ans = numpy.fromstring(d, dtype=np_array.dtype)
    >>> bool((np_array == np_ans).all())
    True

    >>> import ctypes
    >>> typesize = 8
    >>> data = [float(i) for i in range(items)]
    >>> Array = ctypes.c_double * items
    >>> a = Array(*data)
    >>> c = blosc.compress_ptr(ctypes.addressof(a), items, typesize)
    >>> d = blosc.decompress(c)
    >>> import struct
    >>> ans = [struct.unpack('d', d[i:i+typesize])[0] \
            for i in range(0, items*typesize, typesize)]
    >>> data == ans
    True
    """

    _check_address(address)
    if items < 0:
        raise ValueError("items cannot be negative")
    length = items * typesize
    _check_input_length('length', length)
    _check_typesize(typesize)
    _check_shuffle(shuffle)
    _check_clevel(clevel)
    _check_cname(cname)

    return _ext.compress_ptr(address, length, typesize, clevel, shuffle, cname)


def decompress(bytes_like, as_bytearray=False):
    """decompress(bytes_like)

    Decompresses a bytes-like compressed object.

    Parameters
    ----------
    bytes_like : bytes-like object
        The data to be decompressed.  Must be a bytes-like object
        that supports the Python Buffer Protocol, like bytes, bytearray,
        memoryview, or numpy.ndarray.
    as_bytearray : bool, optional
        If this flag is True then the return type will be a bytearray object
        instead of a bytesobject.

    Returns
    -------
    out : str / bytes or bytearray
        The decompressed data in form of a Python str / bytes object.
        If as_bytearray is True then this will be a bytearray object, otherwise
        this will be a str/ bytes object.

    Raises
    ------
    TypeError
        If bytes_like does not support Buffer Protocol

    Examples
    --------

    >>> import array, sys
    >>> a = array.array('i', range(1000*1000))
    >>> a_bytesobj = a.tobytes()
    >>> c_bytesobj = blosc.compress(a_bytesobj, typesize=4)
    >>> a_bytesobj2 = blosc.decompress(c_bytesobj)
    >>> a_bytesobj == a_bytesobj2
    True
    >>> b"" == blosc.decompress(blosc.compress(b"", 1))
    True
    >>> b"1"*7 == blosc.decompress(blosc.compress(b"1"*7, 8))
    True
    >>> type(blosc.decompress(blosc.compress(b"1"*7, 8),
    ...                                      as_bytearray=True)) is bytearray
    True

    """

    return _ext.decompress(bytes_like, as_bytearray)


def decompress_ptr(bytes_like, address):
    """decompress_ptr(bytes_like, address)

    Decompresses a bytes_like compressed object into the memory at address.

    Parameters
    ----------
    bytes_like : bytes-like object
        The data to be decompressed. Must be a bytes-like object
        that supports the Python Buffer Protocol, like bytes, bytearray,
        memoryview, or numpy.ndarray.
    address : int or long
        The address at which to write the decompressed data

    Returns
    -------
    nbytes : int
        the number of bytes written to the buffer

    Raises
    ------
    TypeError
        If bytesobj is not of type bytes or string.
        If address is not of type int or long.

    Notes
    -----
    This function can be used anywhere that a memory address is available in
    Python. For example the Numpy "__array_interface__['data'][0]" construct,
    or when using the ctypes modules.

    Importantly, the user is responsible for making sure that the memory
    address is valid and that the memory pointed to is contiguous and can be
    written to. Passing a non-valid address has a high likelihood of crashing
    the interpreter by segfault.

    Examples
    --------

    >>> import numpy
    >>> items = 7
    >>> np_array = numpy.arange(items)
    >>> c = blosc.compress_ptr(np_array.__array_interface__['data'][0], \
        items, np_array.dtype.itemsize)
    >>> np_ans = numpy.empty(items, dtype=np_array.dtype)
    >>> nbytes = blosc.decompress_ptr(c, np_ans.__array_interface__['data'][0])
    >>> bool((np_array == np_ans).all())
    True
    >>> nbytes == items * np_array.dtype.itemsize
    True

    >>> import ctypes
    >>> typesize = 8
    >>> data = [float(i) for i in range(items)]
    >>> Array = ctypes.c_double * items
    >>> in_array = Array(*data)
    >>> c = blosc.compress_ptr(ctypes.addressof(in_array), items, typesize)
    >>> out_array = ctypes.create_string_buffer(items*typesize)
    >>> nbytes = blosc.decompress_ptr(c, ctypes.addressof(out_array))
    >>> import struct
    >>> ans = [struct.unpack('d', out_array[i:i+typesize])[0] \
            for i in range(0, items*typesize, typesize)]
    >>> data == ans
    True
    >>> nbytes == items * typesize
    True

    """

    _check_byteslike(bytes_like)
    _check_address(address)

    return _ext.decompress_ptr(bytes_like, address)


def pack_array(array, clevel=9, shuffle=blosc.SHUFFLE, cname='blosclz'):
    """pack_array(array[, clevel=9, shuffle=blosc.SHUFFLE, cname='blosclz']])

    Pack (compress) a NumPy array.

    Parameters
    ----------
    array : ndarray
        The NumPy array to be packed.
    clevel : int (optional)
        The compression level from 0 (no compression) to 9
        (maximum compression).  The default is 9.
    shuffle : int (optional)
        The shuffle filter to be activated.  Allowed values are
        blosc.NOSHUFFLE, blosc.SHUFFLE and blosc.BITSHUFFLE.  The
        default is blosc.SHUFFLE.
    cname : string (optional)
        The name of the compressor used internally in Blosc. It can be
        any of the supported by Blosc ('blosclz', 'lz4', 'lz4hc',
        'snappy', 'zlib', 'zstd' and maybe others too). The default is
        'blosclz'.

    Returns
    -------
    out : str / bytes
        The packed array in form of a Python str / bytes object.

    Raises
    ------
    TypeError
        If array does not quack like a numpy ndarray.
    ValueError
        If array.itemsize * array.size is larger than the maximum allowed buffer size.
        If typesize is not within the allowed range.
        If clevel is not within the allowed range.
        If cname is not within the supported compressors.

    Examples
    --------

    >>> import numpy
    >>> a = numpy.arange(1e6)
    >>> parray = blosc.pack_array(a)
    >>> len(parray) < a.size*a.itemsize
    True

    """

    if not (hasattr(array, 'dtype') and hasattr(array, 'shape')):
        # This does not quack like an ndarray
        raise TypeError(
            "only NumPy ndarray objects supported as input")
    itemsize = array.itemsize
    _check_input_length('array size', array.size*itemsize)
    _check_typesize(array.itemsize)
    _check_shuffle(shuffle)
    _check_clevel(clevel)
    _check_cname(cname)

    # Use the fastest pickle available
    pickled_array = pickle.dumps(array, pickle.HIGHEST_PROTOCOL)
    # ... and compress the pickle
    packed_array = compress(pickled_array, itemsize, clevel, shuffle, cname)

    return packed_array


def unpack_array(packed_array, **kwargs):
    """unpack_array(packed_array)

    Unpack (decompress) a packed NumPy array.

    Parameters
    ----------
    packed_array : str / bytes
        The packed array to be decompressed.

    **kwargs : fix_imports / encoding / errors
        Optional parameters that can be passed to the pickle.loads API
        https://docs.python.org/3/library/pickle.html#pickle.loads

    Returns
    -------
    out : ndarray
        The decompressed data in form of a NumPy array.

    Raises
    ------
    TypeError
        If packed_array is not of type bytes or string.

    Examples
    --------

    >>> import numpy
    >>> a = numpy.arange(1e6)
    >>> parray = blosc.pack_array(a)
    >>> len(parray) < a.size*a.itemsize
    True
    >>> a2 = blosc.unpack_array(parray)
    >>> bool(numpy.all(a == a2))
    True
    >>> a = numpy.array(['å', 'ç', 'ø'])
    >>> parray = blosc.pack_array(a)
    >>> a2 = blosc.unpack_array(parray)
    >>> bool(numpy.all(a == a2))
    True
    """

    _check_bytesobj(packed_array)

    # First decompress the pickle
    pickled_array = _ext.decompress(packed_array, False)
    # ... and unpickle

    if kwargs:
        array = pickle.loads(pickled_array, **kwargs)
        if all(isinstance(x, bytes) for x in array.tolist()):
            import numpy
            array = numpy.array([x.decode('utf-8') for x in array.tolist()])
    else:
        array = pickle.loads(pickled_array)

    return array


# For the load tests protocol:
# https://docs.python.org/3/library/unittest.html#load-tests-protocol
def load_tests(loader, tests, pattern):
    import doctest
    tests.addTests(doctest.DocTestSuite())
    return tests

def os_release_pretty_name():
    for p in ('/etc/os-release', '/usr/lib/os-release'):
        try:
            f = open(p)
            for line in f:
                name, _, value = line.rstrip().partition('=')
                if name == 'PRETTY_NAME':
                    if len(value) >= 2 and value[0] in '"\'' and value[0] == value[-1]:
                        value = value[1:-1]
                    return value
        except OSError:
            pass
    return None

def print_versions():
    """Print all the versions of software that python-blosc relies on."""
    import platform
    print("-=" * 38)
    print("python-blosc version: %s" % blosc.__version__)
    print("Blosc version: %s" % blosc.blosclib_version)
    print("Compressors available: %s" % blosc.cnames)
    print("Compressor library versions:")
    for clib in sorted(blosc.clib_versions.keys()):
        print("  %s: %s" % (clib, blosc.clib_versions[clib]))
    print("Python version: %s" % sys.version)
    (sysname, nodename, release, version, machine, processor) = platform.uname()
    print("Platform: %s-%s-%s (%s)" % (sysname, release, machine, version))
    if sysname == "Linux":
        distro = os_release_pretty_name()
        if distro:
            print("Linux dist:", distro)
    if not processor:
        processor = "not recognized"
    print("Processor: %s" % processor)
    print("Byte-ordering: %s" % sys.byteorder)
    print("Detected cores: %s" % blosc.ncores)
    print("Number of threads to use by default: %s" % blosc.nthreads)
    print("-=" * 38)


if __name__ == '__main__':
    # test myself
    import doctest
    print_versions()
    nfail, ntests = doctest.testmod()
    if nfail == 0:
        print("All %d tests passed successfully!" % ntests)
