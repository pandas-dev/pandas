"""
Tests for the pandas.io.common functionalities
"""
import mmap
import os

import pytest

from pandas.compat import StringIO, is_platform_windows
import pandas.util._test_decorators as td

import pandas as pd
import pandas.util.testing as tm

import pandas.io.common as icom


class CustomFSPath(object):
    """For testing fspath on unknown objects"""
    def __init__(self, path):
        self.path = path

    def __fspath__(self):
        return self.path


# Functions that consume a string path and return a string or path-like object
path_types = [str, CustomFSPath]

try:
    from pathlib import Path
    path_types.append(Path)
except ImportError:
    pass

try:
    from py.path import local as LocalPath
    path_types.append(LocalPath)
except ImportError:
    pass

HERE = os.path.abspath(os.path.dirname(__file__))


# https://github.com/cython/cython/issues/1720
@pytest.mark.filterwarnings("ignore:can't resolve package:ImportWarning")
class TestCommonIOCapabilities(object):
    data1 = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

    def test_expand_user(self):
        filename = '~/sometest'
        expanded_name = icom._expand_user(filename)

        assert expanded_name != filename
        assert os.path.isabs(expanded_name)
        assert os.path.expanduser(filename) == expanded_name

    def test_expand_user_normal_path(self):
        filename = '/somefolder/sometest'
        expanded_name = icom._expand_user(filename)

        assert expanded_name == filename
        assert os.path.expanduser(filename) == expanded_name

    @td.skip_if_no('pathlib')
    def test_stringify_path_pathlib(self):
        rel_path = icom._stringify_path(Path('.'))
        assert rel_path == '.'
        redundant_path = icom._stringify_path(Path('foo//bar'))
        assert redundant_path == os.path.join('foo', 'bar')

    @td.skip_if_no('py.path')
    def test_stringify_path_localpath(self):
        path = os.path.join('foo', 'bar')
        abs_path = os.path.abspath(path)
        lpath = LocalPath(path)
        assert icom._stringify_path(lpath) == abs_path

    def test_stringify_path_fspath(self):
        p = CustomFSPath('foo/bar.csv')
        result = icom._stringify_path(p)
        assert result == 'foo/bar.csv'

    @pytest.mark.parametrize('extension,expected', [
        ('', None),
        ('.gz', 'gzip'),
        ('.bz2', 'bz2'),
        ('.zip', 'zip'),
        ('.xz', 'xz'),
    ])
    @pytest.mark.parametrize('path_type', path_types)
    def test_infer_compression_from_path(self, extension, expected, path_type):
        path = path_type('foo/bar.csv' + extension)
        compression = icom._infer_compression(path, compression='infer')
        assert compression == expected

    def test_get_filepath_or_buffer_with_path(self):
        filename = '~/sometest'
        filepath_or_buffer, _, _, should_close = icom.get_filepath_or_buffer(
            filename)
        assert filepath_or_buffer != filename
        assert os.path.isabs(filepath_or_buffer)
        assert os.path.expanduser(filename) == filepath_or_buffer
        assert not should_close

    def test_get_filepath_or_buffer_with_buffer(self):
        input_buffer = StringIO()
        filepath_or_buffer, _, _, should_close = icom.get_filepath_or_buffer(
            input_buffer)
        assert filepath_or_buffer == input_buffer
        assert not should_close

    def test_iterator(self):
        reader = pd.read_csv(StringIO(self.data1), chunksize=1)
        result = pd.concat(reader, ignore_index=True)
        expected = pd.read_csv(StringIO(self.data1))
        tm.assert_frame_equal(result, expected)

        # GH12153
        it = pd.read_csv(StringIO(self.data1), chunksize=1)
        first = next(it)
        tm.assert_frame_equal(first, expected.iloc[[0]])
        tm.assert_frame_equal(pd.concat(it), expected.iloc[1:])

    @pytest.mark.parametrize('reader, module, error_class, fn_ext', [
        (pd.read_csv, 'os', FileNotFoundError, 'csv'),
        (pd.read_fwf, 'os', FileNotFoundError, 'txt'),
        (pd.read_excel, 'xlrd', FileNotFoundError, 'xlsx'),
        (pd.read_feather, 'feather', Exception, 'feather'),
        (pd.read_hdf, 'tables', FileNotFoundError, 'h5'),
        (pd.read_stata, 'os', FileNotFoundError, 'dta'),
        (pd.read_sas, 'os', FileNotFoundError, 'sas7bdat'),
        (pd.read_json, 'os', ValueError, 'json'),
        (pd.read_msgpack, 'os', ValueError, 'mp'),
        (pd.read_pickle, 'os', FileNotFoundError, 'pickle'),
    ])
    def test_read_non_existant(self, reader, module, error_class, fn_ext):
        pytest.importorskip(module)

        path = os.path.join(HERE, 'data', 'does_not_exist.' + fn_ext)
        msg1 = (r"File (b')?.+does_not_exist\.{}'? does not exist"
                .format(fn_ext))
        msg2 = (r"\[Errno 2\] No such file or directory: '.+does_not_exist"
                r"\.{}'").format(fn_ext)
        msg3 = "Expected object or value"
        msg4 = "path_or_buf needs to be a string file path or file-like"
        msg5 = (r"\[Errno 2\] File .+does_not_exist\.{} does not exist:"
                r" '.+does_not_exist\.{}'").format(fn_ext, fn_ext)
        with pytest.raises(error_class, match=r"({}|{}|{}|{}|{})".format(
                msg1, msg2, msg3, msg4, msg5)):
            reader(path)

    @pytest.mark.parametrize('reader, module, error_class, fn_ext', [
        (pd.read_csv, 'os', FileNotFoundError, 'csv'),
        (pd.read_fwf, 'os', FileNotFoundError, 'txt'),
        (pd.read_excel, 'xlrd', FileNotFoundError, 'xlsx'),
        (pd.read_feather, 'feather', Exception, 'feather'),
        (pd.read_hdf, 'tables', FileNotFoundError, 'h5'),
        (pd.read_stata, 'os', FileNotFoundError, 'dta'),
        (pd.read_sas, 'os', FileNotFoundError, 'sas7bdat'),
        (pd.read_json, 'os', ValueError, 'json'),
        (pd.read_msgpack, 'os', ValueError, 'mp'),
        (pd.read_pickle, 'os', FileNotFoundError, 'pickle'),
    ])
    def test_read_expands_user_home_dir(self, reader, module,
                                        error_class, fn_ext, monkeypatch):
        pytest.importorskip(module)

        path = os.path.join('~', 'does_not_exist.' + fn_ext)
        monkeypatch.setattr(icom, '_expand_user',
                            lambda x: os.path.join('foo', x))

        msg1 = (r"File (b')?.+does_not_exist\.{}'? does not exist"
                .format(fn_ext))
        msg2 = (r"\[Errno 2\] No such file or directory:"
                r" '.+does_not_exist\.{}'").format(fn_ext)
        msg3 = "Unexpected character found when decoding 'false'"
        msg4 = "path_or_buf needs to be a string file path or file-like"
        msg5 = (r"\[Errno 2\] File .+does_not_exist\.{} does not exist:"
                r" '.+does_not_exist\.{}'").format(fn_ext, fn_ext)

        with pytest.raises(error_class, match=r"({}|{}|{}|{}|{})".format(
                msg1, msg2, msg3, msg4, msg5)):
            reader(path)

    def test_read_non_existant_read_table(self):
        path = os.path.join(HERE, 'data', 'does_not_exist.' + 'csv')
        msg1 = r"File b'.+does_not_exist\.csv' does not exist"
        msg2 = (r"\[Errno 2\] File .+does_not_exist\.csv does not exist:"
                r" '.+does_not_exist\.csv'")
        with pytest.raises(FileNotFoundError, match=r"({}|{})".format(
                msg1, msg2)):
            with tm.assert_produces_warning(FutureWarning):
                pd.read_table(path)

    @pytest.mark.parametrize('reader, module, path', [
        (pd.read_csv, 'os', ('io', 'data', 'iris.csv')),
        (pd.read_fwf, 'os', ('io', 'data', 'fixed_width_format.txt')),
        (pd.read_excel, 'xlrd', ('io', 'data', 'test1.xlsx')),
        (pd.read_feather, 'feather', ('io', 'data', 'feather-0_3_1.feather')),
        (pd.read_hdf, 'tables', ('io', 'data', 'legacy_hdf',
                                 'datetimetz_object.h5')),
        (pd.read_stata, 'os', ('io', 'data', 'stata10_115.dta')),
        (pd.read_sas, 'os', ('io', 'sas', 'data', 'test1.sas7bdat')),
        (pd.read_json, 'os', ('io', 'json', 'data', 'tsframe_v012.json')),
        (pd.read_msgpack, 'os', ('io', 'msgpack', 'data', 'frame.mp')),
        (pd.read_pickle, 'os', ('io', 'data', 'categorical_0_14_1.pickle')),
    ])
    def test_read_fspath_all(self, reader, module, path, datapath):
        pytest.importorskip(module)
        path = datapath(*path)

        mypath = CustomFSPath(path)
        result = reader(mypath)
        expected = reader(path)

        if path.endswith('.pickle'):
            # categorical
            tm.assert_categorical_equal(result, expected)
        else:
            tm.assert_frame_equal(result, expected)

    def test_read_fspath_all_read_table(self, datapath):
        path = datapath('io', 'data', 'iris.csv')

        mypath = CustomFSPath(path)
        with tm.assert_produces_warning(FutureWarning):
            result = pd.read_table(mypath)
        with tm.assert_produces_warning(FutureWarning):
            expected = pd.read_table(path)

        if path.endswith('.pickle'):
            # categorical
            tm.assert_categorical_equal(result, expected)
        else:
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('writer_name, writer_kwargs, module', [
        ('to_csv', {}, 'os'),
        ('to_excel', {'engine': 'xlwt'}, 'xlwt'),
        ('to_feather', {}, 'feather'),
        ('to_html', {}, 'os'),
        ('to_json', {}, 'os'),
        ('to_latex', {}, 'os'),
        ('to_msgpack', {}, 'os'),
        ('to_pickle', {}, 'os'),
        ('to_stata', {'time_stamp': pd.to_datetime('2019-01-01 00:00')}, 'os'),
    ])
    def test_write_fspath_all(self, writer_name, writer_kwargs, module):
        p1 = tm.ensure_clean('string')
        p2 = tm.ensure_clean('fspath')
        df = pd.DataFrame({"A": [1, 2]})

        with p1 as string, p2 as fspath:
            pytest.importorskip(module)
            mypath = CustomFSPath(fspath)
            writer = getattr(df, writer_name)

            writer(string, **writer_kwargs)
            with open(string, 'rb') as f:
                expected = f.read()

            writer(mypath, **writer_kwargs)
            with open(fspath, 'rb') as f:
                result = f.read()

            assert result == expected

    def test_write_fspath_hdf5(self):
        # Same test as write_fspath_all, except HDF5 files aren't
        # necessarily byte-for-byte identical for a given dataframe, so we'll
        # have to read and compare equality
        pytest.importorskip('tables')

        df = pd.DataFrame({"A": [1, 2]})
        p1 = tm.ensure_clean('string')
        p2 = tm.ensure_clean('fspath')

        with p1 as string, p2 as fspath:
            mypath = CustomFSPath(fspath)
            df.to_hdf(mypath, key='bar')
            df.to_hdf(string, key='bar')

            result = pd.read_hdf(fspath, key='bar')
            expected = pd.read_hdf(string, key='bar')

        tm.assert_frame_equal(result, expected)


@pytest.fixture
def mmap_file(datapath):
    return datapath('io', 'data', 'test_mmap.csv')


class TestMMapWrapper(object):

    def test_constructor_bad_file(self, mmap_file):
        non_file = StringIO('I am not a file')
        non_file.fileno = lambda: -1

        # the error raised is different on Windows
        if is_platform_windows():
            msg = "The parameter is incorrect"
            err = OSError
        else:
            msg = "[Errno 22]"
            err = mmap.error

        with pytest.raises(err, match=msg):
            icom.MMapWrapper(non_file)

        target = open(mmap_file, 'r')
        target.close()

        msg = "I/O operation on closed file"
        with pytest.raises(ValueError, match=msg):
            icom.MMapWrapper(target)

    def test_get_attr(self, mmap_file):
        with open(mmap_file, 'r') as target:
            wrapper = icom.MMapWrapper(target)

        attrs = dir(wrapper.mmap)
        attrs = [attr for attr in attrs
                 if not attr.startswith('__')]
        attrs.append('__next__')

        for attr in attrs:
            assert hasattr(wrapper, attr)

        assert not hasattr(wrapper, 'foo')

    def test_next(self, mmap_file):
        with open(mmap_file, 'r') as target:
            wrapper = icom.MMapWrapper(target)
            lines = target.readlines()

        for line in lines:
            next_line = next(wrapper)
            assert next_line.strip() == line.strip()

        with pytest.raises(StopIteration, match=r'^$'):
            next(wrapper)

    def test_unknown_engine(self):
        with tm.ensure_clean() as path:
            df = tm.makeDataFrame()
            df.to_csv(path)
            with pytest.raises(ValueError, match='Unknown engine'):
                pd.read_csv(path, engine='pyt')
