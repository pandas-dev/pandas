# -*- coding: utf-8 -*-
"""winpty wrapper tests."""

# Standard library imports
import asyncio
import os
import signal
import time
import sys
import re

# Third party imports
import pytest
from flaky import flaky

# Local imports
from winpty.enums import Backend
from winpty.ptyprocess import PtyProcess, which



@pytest.fixture(scope='module', params=['WinPTY', 'ConPTY'])
def pty_fixture(request):
    backend = request.param
    if os.environ.get('CI_RUNNING', None) == '1':
        if backend == 'ConPTY':
            os.environ['CI'] = '1'
            os.environ['CONPTY_CI'] = '1'
        if backend == 'WinPTY':
            os.environ.pop('CI', None)
            os.environ.pop('CONPTY_CI', None)

    backend = getattr(Backend, backend)
    def _pty_factory(cmd=None, env=None):
        cmd = cmd or 'cmd'
        return PtyProcess.spawn(cmd, env=env, backend=backend)
    _pty_factory.backend = request.param
    return _pty_factory


@flaky(max_runs=40, min_passes=1)
def test_read(pty_fixture):
    pty = pty_fixture()
    loc = os.getcwd()
    data = ''
    tries = 0
    while loc not in data and tries < 10:
        try:
            data += pty.read()
        except EOFError:
            pass
        tries += 1
    assert loc in data
    pty.terminate()
    time.sleep(2)


@flaky(max_runs=40, min_passes=1)
def test_write(pty_fixture):
    pty = pty_fixture()

    text = 'Eggs, ham and spam ünicode'
    pty.write(text)

    data = ''
    tries = 0
    while text not in data and tries < 10:
        try:
            data += pty.read()
        except EOFError:
            pass
        tries += 1
    assert text in data
    pty.terminate()


@pytest.mark.xfail(reason="It fails sometimes due to long strings")
@flaky(max_runs=40, min_passes=1)
def test_isalive(pty_fixture):
    pty = pty_fixture()

    pty.write('echo \"foo\"\r\nexit\r\n')
    data = ''
    while True:
        try:
            print('Stuck')
            data += pty.read()
        except EOFError:
            break

    regex = re.compile(".*foo.*")
    assert regex.findall(data)
    assert not pty.isalive()
    pty.terminate()


@pytest.mark.xfail(reason="It fails sometimes due to long strings")
@flaky(max_runs=40, min_passes=1)
def test_readline(pty_fixture):
    env = os.environ.copy()
    env['foo'] = 'bar'
    pty = pty_fixture(env=env)

    # Ensure that the echo print has its own CRLF
    pty.write('cls\r\n')
    pty.write('echo %foo%\r\n')

    data = ''
    tries = 0
    while 'bar' not in data and tries < 10:
        data = pty.readline()
        tries += 1

    assert 'bar' in data

    pty.terminate()


def test_close(pty_fixture):
    pty = pty_fixture()
    pty.close()
    assert not pty.isalive()


def test_flush(pty_fixture):
    pty = pty_fixture()
    pty.flush()
    pty.terminate()


def test_intr(pty_fixture):
    pty = pty_fixture(cmd=[sys.executable, 'import time; time.sleep(10)'])
    pty.sendintr()
    assert pty.wait() != 0


def test_send_control(pty_fixture):
    pty = pty_fixture(cmd=[sys.executable, 'import time; time.sleep(10)'])
    pty.sendcontrol('d')
    assert pty.wait() != 0


@pytest.mark.skipif(which('cat') is None, reason="Requires cat on the PATH")
def test_send_eof(pty_fixture):
    cat = pty_fixture('cat')
    cat.sendeof()
    assert cat.wait() == 0


def test_isatty(pty_fixture):
    pty = pty_fixture()
    assert pty.isatty()
    pty.terminate()
    assert not pty.isatty()


def test_wait(pty_fixture):
    pty = pty_fixture(cmd=[sys.executable, '--version'])
    assert pty.wait() == 0


def test_exit_status(pty_fixture):
    pty = pty_fixture(cmd=[sys.executable])
    pty.write('import sys;sys.exit(1)\r\n')
    pty.wait()
    assert pty.exitstatus == 1


@pytest.mark.timeout(30)
def test_kill_sigterm(pty_fixture):
    pty = pty_fixture()
    pty.write('echo \"foo\"\r\nsleep 1000\r\n')
    pty.read()
    pty.kill(signal.SIGTERM)

    while True:
        try:
            pty.read()
        except EOFError:
            break

    assert not pty.isalive()
    assert pty.exitstatus == signal.SIGTERM


@pytest.mark.timeout(30)
def test_terminate(pty_fixture):
    pty = pty_fixture()
    pty.write('echo \"foo\"\r\nsleep 1000\r\n')
    pty.read()
    pty.terminate()

    while True:
        try:
            pty.read()
        except EOFError:
            break

    assert not pty.isalive()
    assert pty.closed


@pytest.mark.timeout(30)
def test_terminate_force(pty_fixture):
    pty = pty_fixture()
    pty.write('echo \"foo\"\r\nsleep 1000\r\n')
    pty.read()
    pty.terminate(force=True)

    while True:
        try:
            pty.read()
        except EOFError:
            break

    assert not pty.isalive()
    assert pty.closed


def test_terminate_loop(pty_fixture):
    pty = pty_fixture()
    loop = asyncio.SelectorEventLoop()
    asyncio.set_event_loop(loop)

    def reader():
        try:
            data = pty.read()
        except EOFError:
            loop.remove_reader(pty.fd)
            loop.stop()

    loop.add_reader(pty.fd, reader)
    loop.call_soon(pty.write, 'echo \"foo\"\r\nsleep 1000\r\n')
    loop.call_soon(pty.terminate, True)

    try:
        loop.run_forever()
    finally:
        loop.close()

    assert not pty.isalive()
    assert pty.closed


def test_getwinsize(pty_fixture):
    pty = pty_fixture()
    assert pty.getwinsize() == (24, 80)
    pty.terminate()


def test_setwinsize(pty_fixture):
    pty = pty_fixture()
    pty.setwinsize(50, 110)
    assert pty.getwinsize() == (50, 110)
    pty.terminate()

    pty = PtyProcess.spawn('cmd', dimensions=(60, 120))
    assert pty.getwinsize() == (60, 120)
    pty.terminate()

