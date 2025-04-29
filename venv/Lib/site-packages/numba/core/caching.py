"""
Caching mechanism for compiled functions.
"""


from abc import ABCMeta, abstractmethod, abstractproperty
import contextlib
import errno
import hashlib
import inspect
import itertools
import os
import pickle
import sys
import tempfile
import uuid
import warnings

from numba.misc.appdirs import AppDirs
import zipfile
from pathlib import Path

import numba
from numba.core.errors import NumbaWarning
from numba.core.base import BaseContext
from numba.core.codegen import CodeLibrary
from numba.core.compiler import CompileResult
from numba.core import config, compiler
from numba.core.serialize import dumps


def _cache_log(msg, *args):
    if config.DEBUG_CACHE:
        msg = msg % args
        print(msg)


class _Cache(metaclass=ABCMeta):

    @abstractproperty
    def cache_path(self):
        """
        The base filesystem path of this cache (for example its root folder).
        """

    @abstractmethod
    def load_overload(self, sig, target_context):
        """
        Load an overload for the given signature using the target context.
        The saved object must be returned if successful, None if not found
        in the cache.
        """

    @abstractmethod
    def save_overload(self, sig, data):
        """
        Save the overload for the given signature.
        """

    @abstractmethod
    def enable(self):
        """
        Enable the cache.
        """

    @abstractmethod
    def disable(self):
        """
        Disable the cache.
        """

    @abstractmethod
    def flush(self):
        """
        Flush the cache.
        """


class NullCache(_Cache):
    @property
    def cache_path(self):
        return None

    def load_overload(self, sig, target_context):
        pass

    def save_overload(self, sig, cres):
        pass

    def enable(self):
        pass

    def disable(self):
        pass

    def flush(self):
        pass


class _CacheLocator(metaclass=ABCMeta):
    """
    A filesystem locator for caching a given function.
    """

    def ensure_cache_path(self):
        path = self.get_cache_path()
        os.makedirs(path, exist_ok=True)
        # Ensure the directory is writable by trying to write a temporary file
        tempfile.TemporaryFile(dir=path).close()

    @abstractmethod
    def get_cache_path(self):
        """
        Return the directory the function is cached in.
        """

    @abstractmethod
    def get_source_stamp(self):
        """
        Get a timestamp representing the source code's freshness.
        Can return any picklable Python object.
        """

    @abstractmethod
    def get_disambiguator(self):
        """
        Get a string disambiguator for this locator's function.
        It should allow disambiguating different but similarly-named functions.
        """

    @classmethod
    def from_function(cls, py_func, py_file):
        """
        Create a locator instance for the given function located in the
        given file.
        """
        raise NotImplementedError

    @classmethod
    def get_suitable_cache_subpath(cls, py_file):
        """Given the Python file path, compute a suitable path inside the
        cache directory.

        This will reduce a file path that is too long, which can be a problem
        on some operating system (i.e. Windows 7).
        """
        path = os.path.abspath(py_file)
        subpath = os.path.dirname(path)
        parentdir = os.path.split(subpath)[-1]
        # Use SHA1 to reduce path length.
        # Note: windows doesn't like long path.
        hashed = hashlib.sha1(subpath.encode()).hexdigest()
        # Retain parent directory name for easier debugging
        return '_'.join([parentdir, hashed])


class _SourceFileBackedLocatorMixin(object):
    """
    A cache locator mixin for functions which are backed by a well-known
    Python source file.
    """

    def get_source_stamp(self):
        if getattr(sys, 'frozen', False):
            st = os.stat(sys.executable)
        else:
            st = os.stat(self._py_file)
        # We use both timestamp and size as some filesystems only have second
        # granularity.
        return st.st_mtime, st.st_size

    def get_disambiguator(self):
        return str(self._lineno)

    @classmethod
    def from_function(cls, py_func, py_file):
        if not os.path.exists(py_file):
            # Perhaps a placeholder (e.g. "<ipython-XXX>")
            return
        self = cls(py_func, py_file)
        try:
            self.ensure_cache_path()
        except OSError:
            # Cannot ensure the cache directory exists or is writable
            return
        return self


class _UserProvidedCacheLocator(_SourceFileBackedLocatorMixin, _CacheLocator):
    """
    A locator that always point to the user provided directory in
    `numba.config.CACHE_DIR`
    """
    def __init__(self, py_func, py_file):
        self._py_file = py_file
        self._lineno = py_func.__code__.co_firstlineno
        cache_subpath = self.get_suitable_cache_subpath(py_file)
        self._cache_path = os.path.join(config.CACHE_DIR, cache_subpath)

    def get_cache_path(self):
        return self._cache_path

    @classmethod
    def from_function(cls, py_func, py_file):
        if not config.CACHE_DIR:
            return
        parent = super(_UserProvidedCacheLocator, cls)
        return parent.from_function(py_func, py_file)


class _InTreeCacheLocator(_SourceFileBackedLocatorMixin, _CacheLocator):
    """
    A locator for functions backed by a regular Python module with a
    writable __pycache__ directory.
    """

    def __init__(self, py_func, py_file):
        self._py_file = py_file
        self._lineno = py_func.__code__.co_firstlineno
        self._cache_path = os.path.join(os.path.dirname(self._py_file), '__pycache__')

    def get_cache_path(self):
        return self._cache_path


class _UserWideCacheLocator(_SourceFileBackedLocatorMixin, _CacheLocator):
    """
    A locator for functions backed by a regular Python module or a
    frozen executable, cached into a user-wide cache directory.
    """

    def __init__(self, py_func, py_file):
        self._py_file = py_file
        self._lineno = py_func.__code__.co_firstlineno
        appdirs = AppDirs(appname="numba", appauthor=False)
        cache_dir = appdirs.user_cache_dir
        cache_subpath = self.get_suitable_cache_subpath(py_file)
        self._cache_path = os.path.join(cache_dir, cache_subpath)

    def get_cache_path(self):
        return self._cache_path

    @classmethod
    def from_function(cls, py_func, py_file):
        if not (os.path.exists(py_file) or getattr(sys, 'frozen', False)):
            # Perhaps a placeholder (e.g. "<ipython-XXX>")
            # stop function exit if frozen, since it uses a temp placeholder
            return
        self = cls(py_func, py_file)
        try:
            self.ensure_cache_path()
        except OSError:
            # Cannot ensure the cache directory exists or is writable
            return
        return self


class _IPythonCacheLocator(_CacheLocator):
    """
    A locator for functions entered at the IPython prompt (notebook or other).
    """

    def __init__(self, py_func, py_file):
        self._py_file = py_file
        # Note IPython enhances the linecache module to be able to
        # inspect source code of functions defined on the interactive prompt.
        source = inspect.getsource(py_func)
        if isinstance(source, bytes):
            self._bytes_source = source
        else:
            self._bytes_source = source.encode('utf-8')

    def get_cache_path(self):
        # We could also use jupyter_core.paths.jupyter_runtime_dir()
        # In both cases this is a user-wide directory, so we need to
        # be careful when disambiguating if we don't want too many
        # conflicts (see below).
        try:
            from IPython.paths import get_ipython_cache_dir
        except ImportError:
            # older IPython version
            from IPython.utils.path import get_ipython_cache_dir
        return os.path.join(get_ipython_cache_dir(), 'numba_cache')

    def get_source_stamp(self):
        return hashlib.sha256(self._bytes_source).hexdigest()

    def get_disambiguator(self):
        # Heuristic: we don't want too many variants being saved, but
        # we don't want similar named functions (e.g. "f") to compete
        # for the cache, so we hash the first two lines of the function
        # source (usually this will be the @jit decorator + the function
        # signature).
        firstlines = b''.join(self._bytes_source.splitlines(True)[:2])
        return hashlib.sha256(firstlines).hexdigest()[:10]

    @classmethod
    def from_function(cls, py_func, py_file):
        if not (
            py_file.startswith("<ipython-")
            or os.path.basename(os.path.dirname(py_file)).startswith("ipykernel_")
        ):
            return
        self = cls(py_func, py_file)
        try:
            self.ensure_cache_path()
        except OSError:
            # Cannot ensure the cache directory exists
            return
        return self


class _ZipCacheLocator(_SourceFileBackedLocatorMixin, _CacheLocator):
    """
    A locator for functions backed by Python modules within a zip archive.
    """

    def __init__(self, py_func, py_file):
        self._py_file = py_file
        self._lineno = py_func.__code__.co_firstlineno
        self._zip_path, self._internal_path = self._split_zip_path(py_file)
        # We use AppDirs at the moment. A more advanced version of this could also allow
        # a provided `cache_dir`, though that starts to create (cache location x source
        # type) number of cache classes.
        appdirs = AppDirs(appname="numba", appauthor=False)
        cache_dir = appdirs.user_cache_dir
        cache_subpath = self.get_suitable_cache_subpath(py_file)
        self._cache_path = os.path.join(cache_dir, cache_subpath)

    @staticmethod
    def _split_zip_path(py_file):
        path = Path(py_file)
        for i, part in enumerate(path.parts):
            if part.endswith(".zip"):
                zip_path = str(Path(*path.parts[: i + 1]))
                internal_path = str(Path(*path.parts[i + 1 :]))
                return zip_path, internal_path
        raise ValueError("No zip file found in path")

    def get_cache_path(self):
        return self._cache_path

    def get_source_stamp(self):
        st = os.stat(self._zip_path)
        return st.st_mtime, st.st_size

    @classmethod
    def from_function(cls, py_func, py_file):
        if ".zip" not in py_file:
            return None
        return cls(py_func, py_file)

class CacheImpl(metaclass=ABCMeta):
    """
    Provides the core machinery for caching.
    - implement how to serialize and deserialize the data in the cache.
    - control the filename of the cache.
    - provide the cache locator
    """

    _locator_classes = [
        _UserProvidedCacheLocator,
        _InTreeCacheLocator,
        _UserWideCacheLocator,
        _IPythonCacheLocator,
        _ZipCacheLocator,
    ]

    def __init__(self, py_func):
        self._lineno = py_func.__code__.co_firstlineno
        # Get qualname
        try:
            qualname = py_func.__qualname__
        except AttributeError:
            qualname = py_func.__name__
        # Find a locator
        source_path = inspect.getfile(py_func)
        for cls in self._locator_classes:
            locator = cls.from_function(py_func, source_path)
            if locator is not None:
                break
        else:
            raise RuntimeError("cannot cache function %r: no locator available "
                               "for file %r" % (qualname, source_path))
        self._locator = locator
        # Use filename base name as module name to avoid conflict between
        # foo/__init__.py and foo/foo.py
        filename = inspect.getfile(py_func)
        modname = os.path.splitext(os.path.basename(filename))[0]
        fullname = "%s.%s" % (modname, qualname)
        abiflags = getattr(sys, 'abiflags', '')
        self._filename_base = self.get_filename_base(fullname, abiflags)

    def get_filename_base(self, fullname, abiflags):
        # '<' and '>' can appear in the qualname (e.g. '<locals>') but
        # are forbidden in Windows filenames
        fixed_fullname = fullname.replace('<', '').replace('>', '')
        fmt = '%s-%s.py%d%d%s'
        return fmt % (fixed_fullname, self.locator.get_disambiguator(),
                      sys.version_info[0], sys.version_info[1], abiflags)

    @property
    def filename_base(self):
        return self._filename_base

    @property
    def locator(self):
        return self._locator

    @abstractmethod
    def reduce(self, data):
        "Returns the serialized form the data"
        pass

    @abstractmethod
    def rebuild(self, target_context, reduced_data):
        "Returns the de-serialized form of the *reduced_data*"
        pass

    @abstractmethod
    def check_cachable(self, data):
        "Returns True if the given data is cachable; otherwise, returns False."
        pass


class CompileResultCacheImpl(CacheImpl):
    """
    Implements the logic to cache CompileResult objects.
    """

    def reduce(self, cres):
        """
        Returns a serialized CompileResult
        """
        return cres._reduce()

    def rebuild(self, target_context, payload):
        """
        Returns the unserialized CompileResult
        """
        return compiler.CompileResult._rebuild(target_context, *payload)

    def check_cachable(self, cres):
        """
        Check cachability of the given compile result.
        """
        cannot_cache = None
        if any(not x.can_cache for x in cres.lifted):
            cannot_cache = "as it uses lifted code"
        elif cres.library.has_dynamic_globals:
            cannot_cache = ("as it uses dynamic globals "
                            "(such as ctypes pointers and large global arrays)")
        if cannot_cache:
            msg = ('Cannot cache compiled function "%s" %s'
                   % (cres.fndesc.qualname.split('.')[-1], cannot_cache))
            warnings.warn_explicit(msg, NumbaWarning,
                                   self._locator._py_file, self._lineno)
            return False
        return True


class CodeLibraryCacheImpl(CacheImpl):
    """
    Implements the logic to cache CodeLibrary objects.
    """

    _filename_prefix = None  # must be overridden

    def reduce(self, codelib):
        """
        Returns a serialized CodeLibrary
        """
        return codelib.serialize_using_object_code()

    def rebuild(self, target_context, payload):
        """
        Returns the unserialized CodeLibrary
        """
        return target_context.codegen().unserialize_library(payload)

    def check_cachable(self, codelib):
        """
        Check cachability of the given CodeLibrary.
        """
        return not codelib.has_dynamic_globals

    def get_filename_base(self, fullname, abiflags):
        parent = super(CodeLibraryCacheImpl, self)
        res = parent.get_filename_base(fullname, abiflags)
        return '-'.join([self._filename_prefix, res])


class IndexDataCacheFile(object):
    """
    Implements the logic for the index file and data file used by a cache.
    """
    def __init__(self, cache_path, filename_base, source_stamp):
        self._cache_path = cache_path
        self._index_name = '%s.nbi' % (filename_base,)
        self._index_path = os.path.join(self._cache_path, self._index_name)
        self._data_name_pattern = '%s.{number:d}.nbc' % (filename_base,)
        self._source_stamp = source_stamp
        self._version = numba.__version__

    def flush(self):
        self._save_index({})

    def save(self, key, data):
        """
        Save a new cache entry with *key* and *data*.
        """
        overloads = self._load_index()
        try:
            # If key already exists, we will overwrite the file
            data_name = overloads[key]
        except KeyError:
            # Find an available name for the data file
            existing = set(overloads.values())
            for i in itertools.count(1):
                data_name = self._data_name(i)
                if data_name not in existing:
                    break
            overloads[key] = data_name
            self._save_index(overloads)
        self._save_data(data_name, data)

    def load(self, key):
        """
        Load a cache entry with *key*.
        """
        overloads = self._load_index()
        data_name = overloads.get(key)
        if data_name is None:
            return
        try:
            return self._load_data(data_name)
        except OSError:
            # File could have been removed while the index still refers it.
            return

    def _load_index(self):
        """
        Load the cache index and return it as a dictionary (possibly
        empty if cache is empty or obsolete).
        """
        try:
            with open(self._index_path, "rb") as f:
                version = pickle.load(f)
                data = f.read()
        except FileNotFoundError:
            # Index doesn't exist yet?
            return {}
        if version != self._version:
            # This is another version.  Avoid trying to unpickling the
            # rest of the stream, as that may fail.
            return {}
        stamp, overloads = pickle.loads(data)
        _cache_log("[cache] index loaded from %r", self._index_path)
        if stamp != self._source_stamp:
            # Cache is not fresh.  Stale data files will be eventually
            # overwritten, since they are numbered in incrementing order.
            return {}
        else:
            return overloads

    def _save_index(self, overloads):
        data = self._source_stamp, overloads
        data = self._dump(data)
        with self._open_for_write(self._index_path) as f:
            pickle.dump(self._version, f, protocol=-1)
            f.write(data)
        _cache_log("[cache] index saved to %r", self._index_path)

    def _load_data(self, name):
        path = self._data_path(name)
        with open(path, "rb") as f:
            data = f.read()
        tup = pickle.loads(data)
        _cache_log("[cache] data loaded from %r", path)
        return tup

    def _save_data(self, name, data):
        data = self._dump(data)
        path = self._data_path(name)
        with self._open_for_write(path) as f:
            f.write(data)
        _cache_log("[cache] data saved to %r", path)

    def _data_name(self, number):
        return self._data_name_pattern.format(number=number)

    def _data_path(self, name):
        return os.path.join(self._cache_path, name)

    def _dump(self, obj):
        return dumps(obj)

    @contextlib.contextmanager
    def _open_for_write(self, filepath):
        """
        Open *filepath* for writing in a race condition-free way (hopefully).
        uuid4 is used to try and avoid name collisions on a shared filesystem.
        """
        uid = uuid.uuid4().hex[:16]  # avoid long paths
        tmpname = '%s.tmp.%s' % (filepath, uid)
        try:
            with open(tmpname, "wb") as f:
                yield f
            os.replace(tmpname, filepath)
        except Exception:
            # In case of error, remove dangling tmp file
            try:
                os.unlink(tmpname)
            except OSError:
                pass
            raise


class Cache(_Cache):
    """
    A per-function compilation cache.  The cache saves data in separate
    data files and maintains information in an index file.

    There is one index file per function and Python version
    ("function_name-<lineno>.pyXY.nbi") which contains a mapping of
    signatures and architectures to data files.
    It is prefixed by a versioning key and a timestamp of the Python source
    file containing the function.

    There is one data file ("function_name-<lineno>.pyXY.<number>.nbc")
    per function, function signature, target architecture and Python version.

    Separate index and data files per Python version avoid pickle
    compatibility problems.

    Note:
    This contains the driver logic only.  The core logic is provided
    by a subclass of ``CacheImpl`` specified as *_impl_class* in the subclass.
    """

    # The following class variables must be overridden by subclass.
    _impl_class = None

    def __init__(self, py_func):
        self._name = repr(py_func)
        self._py_func = py_func
        self._impl = self._impl_class(py_func)
        self._cache_path = self._impl.locator.get_cache_path()
        # This may be a bit strict but avoids us maintaining a magic number
        source_stamp = self._impl.locator.get_source_stamp()
        filename_base = self._impl.filename_base
        self._cache_file = IndexDataCacheFile(cache_path=self._cache_path,
                                              filename_base=filename_base,
                                              source_stamp=source_stamp)
        self.enable()

    def __repr__(self):
        return "<%s py_func=%r>" % (self.__class__.__name__, self._name)

    @property
    def cache_path(self):
        return self._cache_path

    def enable(self):
        self._enabled = True

    def disable(self):
        self._enabled = False

    def flush(self):
        self._cache_file.flush()

    def load_overload(self, sig, target_context):
        """
        Load and recreate the cached object for the given signature,
        using the *target_context*.
        """
        # Refresh the context to ensure it is initialized
        target_context.refresh()
        with self._guard_against_spurious_io_errors():
            return self._load_overload(sig, target_context)
        # None returned if the `with` block swallows an exception

    def _load_overload(self, sig, target_context):
        if not self._enabled:
            return
        key = self._index_key(sig, target_context.codegen())
        data = self._cache_file.load(key)
        if data is not None:
            data = self._impl.rebuild(target_context, data)
        return data

    def save_overload(self, sig, data):
        """
        Save the data for the given signature in the cache.
        """
        with self._guard_against_spurious_io_errors():
            self._save_overload(sig, data)

    def _save_overload(self, sig, data):
        if not self._enabled:
            return
        if not self._impl.check_cachable(data):
            return
        self._impl.locator.ensure_cache_path()
        key = self._index_key(sig, data.codegen)
        data = self._impl.reduce(data)
        self._cache_file.save(key, data)

    @contextlib.contextmanager
    def _guard_against_spurious_io_errors(self):
        if os.name == 'nt':
            # Guard against permission errors due to accessing the file
            # from several processes (see #2028)
            try:
                yield
            except OSError as e:
                if e.errno != errno.EACCES:
                    raise
        else:
            # No such conditions under non-Windows OSes
            yield

    def _index_key(self, sig, codegen):
        """
        Compute index key for the given signature and codegen.
        It includes a description of the OS, target architecture and hashes of
        the bytecode for the function and, if the function has a __closure__,
        a hash of the cell_contents.
        """
        codebytes = self._py_func.__code__.co_code
        if self._py_func.__closure__ is not None:
            cvars = tuple([x.cell_contents for x in self._py_func.__closure__])
            # Note: cloudpickle serializes a function differently depending
            #       on how the process is launched; e.g. multiprocessing.Process
            cvarbytes = dumps(cvars)
        else:
            cvarbytes = b''

        hasher = lambda x: hashlib.sha256(x).hexdigest()
        return (sig, codegen.magic_tuple(), (hasher(codebytes),
                                             hasher(cvarbytes),))


class FunctionCache(Cache):
    """
    Implements Cache that saves and loads CompileResult objects.
    """
    _impl_class = CompileResultCacheImpl


# Remember used cache filename prefixes.
_lib_cache_prefixes = set([''])


def make_library_cache(prefix):
    """
    Create a Cache class for additional compilation features to cache their
    result for reuse.  The cache is saved in filename pattern like
    in ``FunctionCache`` but with additional *prefix* as specified.
    """
    # avoid cache prefix reuse
    assert prefix not in _lib_cache_prefixes
    _lib_cache_prefixes.add(prefix)

    class CustomCodeLibraryCacheImpl(CodeLibraryCacheImpl):
        _filename_prefix = prefix

    class LibraryCache(Cache):
        """
        Implements Cache that saves and loads CodeLibrary objects for additional
        feature for the specified python function.
        """
        _impl_class = CustomCodeLibraryCacheImpl

    return LibraryCache

