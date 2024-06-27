# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import sys
import shutil
import pickle
import multiprocessing
import threading
import traceback
import datetime

import pytest

from asv import util

WIN = (os.name == 'nt')


def _multiprocessing_raise_processerror(arg):
    try:
        raise util.ProcessError(["a"], 1, "aa", "bb")
    except BaseException as exc:
        # If the following is just 'raise', multiprocessing will hang
        # on Python 2.7.8 due to https://bugs.python.org/issue9400
        raise util.ParallelFailure(str(exc), exc.__class__, traceback.format_exc())


def _multiprocessing_raise_usererror(arg):
    try:
        raise util.UserError("hello")
    except BaseException as exc:
        raise util.ParallelFailure(str(exc), exc.__class__, traceback.format_exc())


@pytest.mark.timeout(30)
def test_parallelfailure():
    # Check the workaround for https://bugs.python.org/issue9400 works

    if WIN and os.path.basename(sys.argv[0]).lower().startswith('py.test'):
        # Multiprocessing in spawn mode can result to problems with py.test
        pytest.skip("Multiprocessing spawn mode on Windows not safe to run "
                    "from py.test runner.")

    # The exception class must be pickleable
    exc = util.ParallelFailure("test", Exception, "something")
    exc2 = pickle.loads(pickle.dumps(exc))
    assert exc.message == exc2.message
    assert exc.exc_cls == exc2.exc_cls
    assert exc.traceback_str == exc2.traceback_str
    assert str(exc) == "Exception: test\n    something"

    # Check multiprocessing does not hang (it would hang on Python
    # 2.7.8 if the 'raise utill.ParallelFailure ...' above is changed
    # to just 'raise')
    pool = multiprocessing.Pool(4)
    try:
        pool.map(_multiprocessing_raise_processerror, range(10))
    except util.ParallelFailure:
        pass
    finally:
        pool.close()

    # Check reraising UserError
    pool = multiprocessing.Pool(4)
    try:
        try:
            pool.map(_multiprocessing_raise_usererror, range(10))
        except util.ParallelFailure as exc:
            exc.reraise()
        finally:
            pool.close()
        assert False
    except util.UserError:
        # OK
        pass


def test_which_path(tmpdir):
    dirname = os.path.abspath(os.path.join(str(tmpdir), 'name with spaces'))
    fn = 'asv_test_exe_1234.exe'
    fn2 = 'asv_test_exe_4321.bat'

    os.makedirs(dirname)
    shutil.copyfile(sys.executable, os.path.join(dirname, fn))
    shutil.copyfile(sys.executable, os.path.join(dirname, fn2))

    old_path = os.environ.get('PATH', '')
    try:
        if WIN:
            os.environ['PATH'] = old_path + os.pathsep + '"' + dirname + '"'
            util.which('asv_test_exe_1234')
            util.which('asv_test_exe_1234.exe')
            util.which('asv_test_exe_4321')
            util.which('asv_test_exe_4321.bat')

        os.environ['PATH'] = old_path + os.pathsep + dirname
        util.which('asv_test_exe_1234.exe')
        util.which('asv_test_exe_4321.bat')
        if WIN:
            util.which('asv_test_exe_1234')
            util.which('asv_test_exe_4321')

        # Check the paths= argument
        util.which('asv_test_exe_1234.exe', paths=[dirname])
        util.which('asv_test_exe_4321.bat', paths=[dirname])

        # Check non-existent files
        with pytest.raises(OSError):
            util.which('nonexistent.exe', paths=[dirname])
    finally:
        os.environ['PATH'] = old_path


def test_write_load_json(tmpdir):
    data = {
        'a': 1,
        'b': 2,
        'c': 3
    }
    orig_data = dict(data)

    filename = os.path.join(str(tmpdir), 'test.json')

    util.write_json(filename, data)
    data2 = util.load_json(filename)
    assert data == orig_data
    assert data2 == orig_data

    util.write_json(filename, data, 3)
    data2 = util.load_json(filename, 3)
    assert data == orig_data
    assert data2 == orig_data

    # Wrong API version must fail to load
    with pytest.raises(util.UserError):
        util.load_json(filename, 2)
    with pytest.raises(util.UserError):
        util.load_json(filename, 4)
    util.write_json(filename, data)
    with pytest.raises(util.UserError):
        util.load_json(filename, 3)


def test_human_float():
    items = [
        # (expected, value, significant, truncate_small, significant_zeros, reference_value)

        # significant
        ("1", 1.2345, 1),
        ("1.2", 1.2345, 2),
        ("1.23", 1.2345, 3),
        ("100", 123.45, 1),
        ("120", 123.45, 2),
        ("123", 123.45, 3),
        ("123.5", 123.45, 4),
        ("0.001", 0.0012346, 1),
        ("0.001235", 0.0012346, 4),

        # significant zeros
        ("0.001", 0.001, 1, None, True),
        ("0.0010", 0.001, 2, None, True),
        ("0.00100", 0.001, 3, None, True),
        ("1", 1, 1, None, True),
        ("1.0", 1, 2, None, True),
        ("1.00", 1, 3, None, True),

        # truncate small
        ("0", 0.001, 2, 0),
        ("0", 0.001, 2, 1),
        ("0.001", 0.001, 2, 2),

        # non-finite
        ("inf", float('inf'), 1),
        ("-inf", -float('inf'), 1),
        ("nan", float('nan'), 1),

        # negative
        ("-1", -1.2345, 1),
        ("-0.00100", -0.001, 3, None, True),
        ("-0", -0.001, 2, 1),
        ("-0.001", -0.001, 2, 2),
    ]

    for item in items:
        expected = item[0]
        got = util.human_float(*item[1:])
        assert got == expected, item


def test_human_time():
    items = [
        # (expected, value, err)

        # scales
        ("1.00ns", 1e-9),
        ("1.10μs", 1.1e-6),
        ("1.12ms", 1.12e-3),
        ("1.12s", 1.123),
        ("1.13s", 1.126),
        ("1.00m", 60),
        ("2.00h", 3600 * 2),
        ("0s", 0),
        ("n/a", float("nan")),

        # err
        ("1.00±1ns", 1e-9, 1e-9),
        ("1.00±0.1ns", 1e-9, 0.1e-9),
        ("1.00±0.01ns", 1e-9, 0.01e-9),
        ("1.00±0.01ns", 1e-9, 0.006e-9),
        ("1.00±0ns", 1e-9, 0.001e-9),
        ("1.00±1000000ns", 1e-9, 1e-3),
        ("0±1s", 0, 1),
        ("0±1ms", 0, 1e-3),
        ("0±0s", 0, 0),
    ]

    for item in items:
        expected = item[0]
        got = util.human_time(*item[1:])
        assert got == expected, item
        got = util.human_value(item[1], 'seconds', *item[2:])
        assert got == expected, item


def test_human_file_size():
    items = [
        # (expected, value, err)

        # scales
        ("1", 1),
        ("999", 999),
        ("1k", 1000),
        ("1.1M", 1.1e6),
        ("1.12G", 1.12e9),
        ("1.12T", 1.123e12),

        # err
        ("1±2", 1, 2),
        ("1±0.1k", 1e3, 123),
        ("12.3±4M", 12.34e6, 4321e3),
    ]

    for item in items:
        expected = item[0]
        got = util.human_file_size(*item[1:])
        assert got == expected, item
        got = util.human_value(item[1], 'bytes', *item[2:])
        assert got == expected, item


def test_parse_human_time():
    items = [
        # (value, expected)
        ("1", 60 * 60 * 24),
        ("1h", 60 * 60),
        ("1w", 60 * 60 * 24 * 7),
    ]

    for value, expected in items:
        result = util.parse_human_time(value)
        assert result == expected

    bad_items = [
        "1:",
        ".",
        "1x",
    ]

    for value in bad_items:
        with pytest.raises(ValueError):
            util.parse_human_time(value)


def test_is_main_thread():
    if sys.version_info[0] >= 3:
        # NB: the test itself doesn't necessarily run in main thread...
        is_main = (threading.current_thread() == threading.main_thread())
        assert util.is_main_thread() == is_main

    results = []

    def worker():
        results.append(util.is_main_thread())

    thread = threading.Thread(target=worker)
    thread.start()
    thread.join()

    assert results == [False]


def test_json_non_ascii(tmpdir):
    non_ascii_data = [{'😼': '難', 'ä': 3}]

    fn = os.path.join(str(tmpdir), "nonascii.json")
    util.write_json(fn, non_ascii_data)
    data = util.load_json(fn)

    assert data == non_ascii_data


def test_interpolate_command():
    good_items = [
        ('python {inputs}', dict(inputs='9'),
         ['python', '9'], {}, {0}, None),

        ('python "{inputs}"', dict(inputs='9'),
         ['python', '9'], {}, {0}, None),

        ('python {inputs}', dict(inputs=''),
         ['python', ''], {}, {0}, None),

        ('HELLO="asd" python "{inputs}"', dict(inputs='9'),
         ['python', '9'], {'HELLO': 'asd'}, {0}, None),

        ('HELLO="asd" return-code=any python "{inputs}"', dict(inputs='9'),
         ['python', '9'], {'HELLO': 'asd'}, None, None),

        ('HELLO="asd" return-code=255 python "{inputs}"', dict(inputs='9'),
         ['python', '9'], {'HELLO': 'asd'}, {255}, None),

        ('HELLO="asd" return-code=255 python "{inputs}"', dict(inputs='9'),
         ['python', '9'], {'HELLO': 'asd'}, {255}, None),

        ('HELLO="asd" in-dir="{somedir}" python', dict(somedir='dir'),
         ['python'], {'HELLO': 'asd'}, {0}, 'dir'),
    ]

    bad_items = [
        ('python {foo}', {}),
        ('HELLO={foo} python', {}),
        ('return-code=none python', {}),
        ('return-code= python', {}),
        ('return-code=, python', {}),
        ('return-code=1,,2 python', {}),
        ('return-code=1 return-code=2 python', {}),
        ('in-dir=a in-dir=b python', {}),
    ]

    for value, variables, e_parts, e_env, e_codes, e_cwd in good_items:
        parts, env, codes, cwd = util.interpolate_command(value, variables)
        assert parts == e_parts
        assert env == e_env
        assert codes == e_codes
        assert cwd == e_cwd

    for value, variables in bad_items:
        with pytest.raises(util.UserError):
            util.interpolate_command(value, variables)


def test_datetime_to_js_timestamp():
    tss = [0, 0.5, -0.5, 12345.6789, -12345.6789,
           1535910708.7767508]
    for ts in tss:
        t = datetime.datetime.fromtimestamp(ts, tz=datetime.timezone.utc)
        ts2 = util.datetime_to_js_timestamp(t)
        assert abs(ts * 1000 - ts2) <= 0.5

    # Check sub-second precision
    ms = 50
    ts = datetime.datetime(
        1970, 1, 1, 0, 0, 0, 1000 * ms,
        tzinfo = datetime.timezone.utc
    )
    assert util.datetime_to_js_timestamp(ts) == ms

    # Check rounding
    ts = datetime.datetime(
        1970, 1, 1, 0, 0, 0, 500,
        tzinfo = datetime.timezone.utc
    )
    assert util.datetime_to_js_timestamp(ts) == 1
    ts = datetime.datetime(
        1970, 1, 1, 0, 0, 0, 499,
        tzinfo = datetime.timezone.utc
    )
    assert util.datetime_to_js_timestamp(ts) == 0


def test_datetime_to_timestamp():
    tss = [0, 0.5, -0.5, 12345.6789, -12345.6789,
           1535910708.7767508]
    for ts in tss:
        t = datetime.datetime.fromtimestamp(ts,
                                            tz=datetime.timezone.utc)
        ts2 = util.datetime_to_timestamp(t)
        assert abs(ts - ts2) <= 0.5

    # Check rounding
    ts = datetime.datetime(
        1970, 1, 1, 0, 0, 0, 500000,
        tzinfo = datetime.timezone.utc
    )
    assert util.datetime_to_timestamp(ts) == 1
    ts = datetime.datetime(
        1970, 1, 1, 0, 0, 0, 500000 - 1,
        tzinfo = datetime.timezone.utc
    )
    assert util.datetime_to_timestamp(ts) == 0


def test_check_output_exit_code(capsys):
    with pytest.raises(util.ProcessError):
        util.check_output([sys.executable, '-c', 'import sys; sys.exit(1)'])
    out, err = capsys.readouterr()
    assert '(exit status 1)' in out


def test_geom_mean_na():
    for x in [[1, 2, -3], [1, 2, 3], [3, 1, 3, None, None]]:
        expected = abs(x[0] * x[1] * x[2])**(1 / 3)
        assert abs(util.geom_mean_na(x) - expected) < 1e-10

def mock_pip_caller(args):
    return args

@pytest.mark.parametrize("declaration, expected", [
    # Basic package name
    ("numpy", {
        'pkgname': 'numpy', 'specification': None,
        'flags': [], 'is_editable': False, 'path': None
    }),

    # Version with '=='
    ("numpy==1.20.0", {
        'pkgname': 'numpy', 'specification': '==1.20.0',
        'flags': [], 'is_editable': False, 'path': None
    }),

    # Other specifiers
    ("numpy>=1.20.0", {
        'pkgname': 'numpy', 'specification': '>=1.20.0',
        'flags': [], 'is_editable': False, 'path': None
    }),

    ("numpy<=1.20.0", {
        'pkgname': 'numpy', 'specification': '<=1.20.0',
        'flags': [], 'is_editable': False, 'path': None
    }),

    # Complex version specifiers
    ("numpy>=1.15,<2.0", {
        'pkgname': 'numpy', 'specification': '>=1.15,<2.0',
        'flags': [], 'is_editable': False, 'path': None
    }),

    # Flags
    ("--no-cache-dir numpy", {
        'pkgname': 'numpy', 'specification': None,
        'flags': ['--no-cache-dir'], 'is_editable': False, 'path': None
    }),

    # Editable installations
    ("-e ./numpy-dir", {
        'pkgname': None, 'specification': None,
        'flags': ['-e'], 'is_editable': True, 'path': './numpy-dir'
    }),

    # Local paths without -e
    ("./numpy-dir", {
        'pkgname': None, 'specification': None,
        'flags': [], 'is_editable': False, 'path': './numpy-dir'
    }),

    # Package with dash
    ("my-package", {
        'pkgname': 'my-package', 'specification': None,
        'flags': [], 'is_editable': False, 'path': None
    }),

    # More complex declarations
    ("--no-cache-dir -e ./numpy-dir", {
        'pkgname': None, 'specification': None,
        'flags': ['--no-cache-dir', '-e'], 'is_editable': True, 'path': './numpy-dir'
    }),

    # Direct git installations
    ("git+https://github.com/user/repo.git", {
        'pkgname': None, 'specification': None,
        'flags': [], 'is_editable': False, 'path': 'git+https://github.com/user/repo.git'
    }),
    # Editable git installations with #egg= suffix
    ("-e git+https://github.com/user/repo.git#egg=repo", {
        'pkgname': 'repo', 'specification': None,
        'flags': ['-e'], 'is_editable': True, 'path': 'git+https://github.com/user/repo.git'
    }),
    # Flags with values
    ("--install-option=\"--prefix=/my/path\" numpy", {
        'pkgname': 'numpy', 'specification': None,
        'flags': ['--install-option=\"--prefix=/my/path\"'], 'is_editable': False, 'path': None
    }),
])
def test_parsed_pip_declaration(declaration, expected):
    parsed = util.ParsedPipDeclaration(declaration)
    for key, value in expected.items():
        assert getattr(parsed, key) == value, \
            f"Expected {key} to be {value}, but got {getattr(parsed, key)}"

# Mock pip_caller
def pip_caller(args):
    return args

@pytest.mark.parametrize(
    "declaration, expected_result",
    [
        # Test with a simple package name
        ("numpy", ["install", "-v", "--upgrade", "numpy"]),

        # Test with a package and version specification
        ("numpy==1.18.5", ["install", "-v", "--upgrade", "numpy==1.18.5"]),

        # Test with a git URL
        ("git+https://github.com/numpy/numpy.git#egg=numpy",
         ["install", "-v", "--upgrade", "git+https://github.com/numpy/numpy.git"]),

        # Test with a local path
        ("./localpackage/", ["install", "-v", "--upgrade", "./localpackage/"]),

        # Test with flags
        ("numpy --install-option=\"--prefix=/my/path\"",
         ["install", "-v", "--upgrade", "--install-option=\"--prefix=/my/path\"", "numpy"]),

        # Test case for multiple version specifiers
        ("numpy>=1.18.5,<=1.20.0", ["install", "-v", "--upgrade", "numpy>=1.18.5,<=1.20.0"]),

        # Test cases for version without specifier
        ("numpy 1.23", ["install", "-v", "--upgrade", "numpy==1.23"]),
        ("fakepy_1 0.14", ["install", "-v", "--upgrade", "fakepy_1==0.14"]),
    ]
)
def test_construct_pip_call(declaration, expected_result):
    parsed_declaration = util.ParsedPipDeclaration(declaration)
    result = util.construct_pip_call(pip_caller, parsed_declaration)
    assert result() == expected_result
