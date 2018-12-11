#!/bin/env python
import os
import subprocess
import sys
import random
import tempfile
import time
import warnings
try:
    from urllib.request import urlretrieve
except ImportError:  # py2
    from urllib import urlretrieve


def set_environ(pattern, locale):
    """
    Set environment variables needed for running the tests.
    """
    # Workaround for pytest-xdist flaky collection order
    # https://github.com/pytest-dev/pytest/issues/920
    # https://github.com/pytest-dev/pytest/issues/1075
    os.environ['PYTHONHASHSEED'] = str(random.randint(1, 4294967295))

    if locale:
        os.environ['LC_ALL'] = os.environ['LANG'] = locale
        import pandas
        pandas_locale = pandas.get_option('display.encoding')
        if pandas_locale != locale:
            # TODO raise exception instead of warning when
            # https://github.com/pandas-dev/pandas/issues/23923 is fixed
            warnings.warn(('pandas could not detect the locale. '
                           'System locale: {}, '
                           'pandas detected: {}').format(locale,
                                                         pandas_locale))

    if 'not network' in pattern:
        os.environ['http_proxy'] = os.environ['https_proxy'] = 'http://1.2.3.4'


def pytest_command(pattern, coverage_file):
    """
    Build and return the pytest command to run.
    """
    cmd = ['pytest', '--junitxml=test-data.xml']

    if pattern:
        cmd += ['-m', pattern]

    if coverage_file:
        cmd += ['--cov=pandas', '--cov-report=xml:{}'.format(coverage_file)]

    test_jobs = os.environ.get('TESTS_JOBS', 0)
    if test_jobs:
        cmd += ['-n', str(test_jobs), '--dist', 'loadfile']

    if os.environ.get('WARNINGS_ARE_ERRORS'):
        cmd += ['-W', 'error']

    return cmd + ['pandas']


def upload_coverage(coverage_file):
    """
    Download codecov.io script and run it to upload coverage for coverage_file.
    """
    script_fname = os.path.join(os.path.dirname(coverage_file),
                                'codecov_script.sh')
    urlretrieve('https://codecov.io/bash', script_fname)
    upload_coverage_cmd = ['bash',
                           script_fname,
                           '-Z',
                           '-c',
                           '-f',
                           coverage_file]
    sys.stderr.write('{}\n'.format(' '.join(upload_coverage_cmd)))
    subprocess.check_call(upload_coverage_cmd.split())
    os.remove(script_fname)
    os.remove(coverage_file)


def run_tests(pattern, locale=None, coverage_file=False):
    """
    Run tests with the specified environment.

    Parameters
    ----------
    pattern : str
        Tests to execute based on pytest markers (e.g. "slow and not network").
    locale : str, optional
        Locale to use instead of the system defaule (e.g. "it_IT.UTF8").
    coverage_file : str, optional
        If provided, the file path where to save the coverage.
    """
    set_environ(pattern, locale)
    pytest_cmd = pytest_command(pattern, coverage_file)
    sys.stderr.write('{}\n'.format(' '.join(pytest_cmd)))
    start = time.time()
    subprocess.check_call(pytest_cmd)
    tests_run_in_seconds = int(time.time() - start)
    sys.stdout.write('Tests run in {} seconds\n'.format(tests_run_in_seconds))

    if coverage_file:
        upload_coverage(coverage_file)


if __name__ == '__main__':
    pattern = os.environ.get('PATTERN', '')
    locale = os.environ.get('LOCALE_OVERRIDE')
    coverage_file = None
    if os.environ.get('COVERAGE', '') != '':
        if sys.platform == 'win32':
            raise RuntimeError('Coverage can not be uploaded from Windows')
        coverage_file = os.path.join(tempfile.gettempdir(),
                                     'pandas-coverage.xml')
    run_tests(pattern, locale, coverage_file)
