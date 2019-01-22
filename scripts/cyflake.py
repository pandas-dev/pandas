import contextlib
import os
import re
import subprocess
import sys
import tempfile


def check_file(path):
    """
    Run a flake8-like check on the cython file at the given path.

    Parameters
    ----------
    path : str

    Returns
    -------
    return_code : int
    """
    with open(path, 'rb') as fd:
        content = fd.read().decode('utf-8')

    py_content = clean_cy_content(content)

    fname = os.path.split(path)[1]

    with ensure_clean(fname) as temp_path:
        with open(temp_path, 'wb') as fd:
            fd.write(py_content.encode('utf-8'))

        rc = call_flake8(temp_path, path)

    return rc


@contextlib.contextmanager
def ensure_clean(filename):
    """
    A poor-man's version of pandas.util.testing.ensure_clean
    """
    fd, filename = tempfile.mkstemp(suffix=filename)

    try:
        yield filename
    finally:
        try:
            os.close(fd)
        except Exception:
            print("Couldn't close file descriptor: {fdesc} (file: {fname})"
                  .format(fdesc=fd, fname=filename))
        try:
            if os.path.exists(filename):
                os.remove(filename)
        except Exception as e:
            print("Exception on removing file: {error}".format(error=e))


def clean_cy_content(content):
    """
    For cython code that we cannot make into valid python, try to edit it
    into something that the linter will recognize.

    Parameters
    ----------
    content : unicode

    Returns
    -------
    unicode
    """

    # Note: this may mess up subsequent lines indentation alignment if there
    #  are multi-line cimports
    content = re.sub(u'^cimport ', u'import ', content)
    content = re.sub(r'(?<=\s)cimport ', u'import ', content)
    return content


def call_flake8(temp_path, orig_path):
    """
    Wrapper to call flake8 on the file at the given temp_path, editing any
    messages to point to the original file's path.

    Parameters
    ----------
    temp_path : str
    orig_path : str

    Returns
    -------
    return_code : int
    """
    p = subprocess.Popen(['flake8', temp_path],
                         close_fds=True,
                         stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    stdout, stderr = p.communicate()

    # Edit the messages to include the original path
    temp_path = temp_path.encode('utf-8')
    orig_path = orig_path.encode('utf-8')
    stdout = stdout.replace(temp_path, orig_path)
    stderr = stderr.replace(temp_path, orig_path)

    # TODO: better to just print?
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

    return p.returncode


if __name__ == "__main__":
    args = sys.argv[1:]
    rc = 0
    for path in args:
        # no validation, we just assume these are paths
        rc2 = check_file(path)
        rc = rc or rc2

    sys.exit(rc)
