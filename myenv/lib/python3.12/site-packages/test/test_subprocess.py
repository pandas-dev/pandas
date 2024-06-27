# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys
import time

import pytest

from asv import util


@pytest.mark.flaky(reruns=3, reruns_delay=5)
def test_timeout():
    timeout_codes = []
    timeout_codes.append(r"""
import sys
import time

sys.stdout.write("Stdout before waiting\n")
sys.stderr.write("Stderr before waiting\n")
sys.stdout.flush()
sys.stderr.flush()
time.sleep(60)
sys.stdout.write("Stdout after waiting\n")
sys.stderr.write("Stderr after waiting\n")
    """)

    # Another example, where timeout is due to a hanging sub-subprocess
    timeout_codes.append(r"""
import sys
import subprocess

sys.stdout.write("Stdout before waiting\n")
sys.stderr.write("Stderr before waiting\n")
sys.stdout.flush()
sys.stderr.flush()
subprocess.call([
    sys.executable,
    "-c",
    "import sys, subprocess; subprocess.call("
    "[sys.executable, '-c', 'import time; time.sleep(360)'])"
])
sys.stdout.write("Stdout after waiting\n")
sys.stderr.write("Stderr after waiting\n")
    """)

    for timeout_code in timeout_codes:
        t = time.time()
        try:
            util.check_output([
                sys.executable, "-c", timeout_code], timeout=1)
        except util.ProcessError as e:
            assert len(e.stdout.strip().split('\n')) == 1
            assert len(e.stderr.strip().split('\n')) == 1
            print(e.stdout)
            assert e.stdout.strip() == "Stdout before waiting"
            assert e.stderr.strip() == "Stderr before waiting"
            assert e.retcode == util.TIMEOUT_RETCODE
            assert "timed out" in str(e)
        else:
            assert False, "Expected timeout exception"
        # Make sure the timeout is triggered in a sufficiently short amount of time
        assert time.time() - t < 5.0


def test_exception():
    code = r"""
import sys
sys.stdout.write("Stdout before error\n")
sys.stderr.write("Stderr before error\n")
sys.exit(1)
"""
    try:
        util.check_output([
            sys.executable, "-c", code])
    except util.ProcessError as e:
        assert len(e.stdout.strip().split('\n')) == 1
        err = [x for x in e.stderr.strip().split('\n')
               if not x.startswith('Coverage')]
        assert len(err) == 1
        assert e.stdout.strip() == "Stdout before error"
        assert err[0] == "Stderr before error"
        assert e.retcode == 1
        assert "returned non-zero exit status 1" in str(e)
    else:
        assert False, "Expected exception"


@pytest.mark.flaky(reruns=3, reruns_delay=5)
def test_output_timeout():
    # Check that timeout is determined based on last output, not based
    # on start time.
    code = r"""
import time
import sys
for j in range(3):
    sys.stdout.write('.')
    sys.stdout.flush()
    time.sleep(1.0)
"""
    output = util.check_output([sys.executable, "-c", code], timeout=1.5)
    assert output == '.' * 3

    try:
        util.check_output([sys.executable, "-c", code], timeout=0.5)
    except util.ProcessError as e:
        assert e.retcode == util.TIMEOUT_RETCODE
    else:
        assert False, "Expected exception"


def test_env():
    code = r"""
import os
print(os.environ['TEST_ASV_FOO'])
print(os.environ['TEST_ASV_BAR'])
"""
    env = os.environ.copy()
    env['TEST_ASV_FOO'] = 'foo'
    # Force unicode string on Python 2
    env['TEST_ASV_BAR'] = u'bar'
    output = util.check_output([sys.executable, "-c", code], env=env)
    assert output.splitlines() == ['foo', 'bar']


def test_no_timeout():
    # Check that timeout=None is allowed.
    code = "import time; time.sleep(0.05)"
    out, err, retcode = util.check_output([sys.executable, "-c", code], timeout=None,
                                          return_stderr=True)
    assert out == ''
    assert err == ''
    assert retcode == 0


def test_stderr_redirect():
    # Check redirecting stderr to stdout works
    code = ("import sys;"
            "sys.stdout.write('OUT\\n');"
            "sys.stdout.flush();"
            "sys.stderr.write('ERR\\n')")
    out = util.check_output([sys.executable, "-c", code], redirect_stderr=True)
    assert out.splitlines() == ['OUT', 'ERR']
    out, err, retcode = util.check_output([sys.executable, "-c", code],
                                          return_stderr=True, redirect_stderr=True)
    assert out.splitlines() == ['OUT', 'ERR']
    assert err == ''
    assert retcode == 0


def test_popen():
    # Check that timeout=None is allowed.
    popen = util.check_output([sys.executable, "-c", "pass"], return_popen=True)
    popen.wait()

    assert popen.returncode == 0

    # close handles
    popen.communicate()


def test_large_output():
    # More data than a pipe buffer can hold
    data = util.check_output([sys.executable, "-c",
                              "import sys; [sys.stdout.write('x'*1000) for j in range(5000)]"])
    assert data == 'x' * 5000000


# This *does* seem to work, only seems untestable somehow...
# def test_dots(capsys):
#     code = r"""
# import sys
# import time
# for i in range(100):
#     sys.stdout.write("Line {0}\n".format(i))
#     sys.stdout.flush()
#     time.sleep(0.001)
# """
#     util.check_output([sys.executable, "-c", code])

#     out, err = capsys.readouterr()

#     assert out == '.' * 100
