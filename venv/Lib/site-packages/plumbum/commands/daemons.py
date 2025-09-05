from __future__ import annotations

import contextlib
import errno
import os
import signal
import subprocess
import sys
import time
import traceback

from plumbum.commands.processes import ProcessExecutionError


class _fake_lock:
    """Needed to allow normal os.exit() to work without error"""

    @staticmethod
    def acquire(_):
        return True

    @staticmethod
    def release():
        pass


def posix_daemonize(command, cwd, stdout=None, stderr=None, append=True):
    if stdout is None:
        stdout = os.devnull
    if stderr is None:
        stderr = stdout

    MAX_SIZE = 16384
    rfd, wfd = os.pipe()
    argv = command.formulate()
    firstpid = os.fork()
    if firstpid == 0:
        # first child: become session leader
        os.close(rfd)
        rc = 0
        try:
            os.setsid()
            os.umask(0)
            stdin = open(os.devnull, encoding="utf-8")
            stdout = open(stdout, "a" if append else "w", encoding="utf-8")
            stderr = open(stderr, "a" if append else "w", encoding="utf-8")
            signal.signal(signal.SIGHUP, signal.SIG_IGN)
            proc = command.popen(
                cwd=cwd,
                close_fds=True,
                stdin=stdin.fileno(),
                stdout=stdout.fileno(),
                stderr=stderr.fileno(),
            )
            os.write(wfd, str(proc.pid).encode("utf8"))
        except Exception:
            rc = 1
            tbtext = "".join(traceback.format_exception(*sys.exc_info()))[-MAX_SIZE:]
            os.write(wfd, tbtext.encode("utf8"))
        finally:
            os.close(wfd)
            os._exit(rc)

    # wait for first child to die
    os.close(wfd)
    _, rc = os.waitpid(firstpid, 0)
    output = os.read(rfd, MAX_SIZE)
    os.close(rfd)
    with contextlib.suppress(UnicodeError):
        output = output.decode("utf8")
    if rc == 0 and output.isdigit():
        secondpid = int(output)
    else:
        raise ProcessExecutionError(argv, rc, "", output)
    proc = subprocess.Popen.__new__(subprocess.Popen)
    proc._child_created = True
    proc.returncode = None
    proc.stdout = None
    proc.stdin = None
    proc.stderr = None
    proc.pid = secondpid
    proc.universal_newlines = False
    proc._input = None
    proc._waitpid_lock = _fake_lock()
    proc._communication_started = False
    proc.args = argv
    proc.argv = argv

    def poll(self=proc):
        if self.returncode is None:
            try:
                os.kill(self.pid, 0)
            except OSError as ex:
                if ex.errno == errno.ESRCH:
                    # process does not exist
                    self.returncode = 0
                else:
                    raise
        return self.returncode

    def wait(self=proc):
        while self.returncode is None:
            if self.poll() is None:
                time.sleep(0.5)
        return proc.returncode

    proc.poll = poll
    proc.wait = wait
    return proc


def win32_daemonize(command, cwd, stdout=None, stderr=None, append=True):
    if stdout is None:
        stdout = os.devnull
    if stderr is None:
        stderr = stdout
    DETACHED_PROCESS = 0x00000008
    stdin = open(os.devnull, encoding="utf-8")
    stdout = open(stdout, "a" if append else "w", encoding="utf-8")
    stderr = open(stderr, "a" if append else "w", encoding="utf-8")
    return command.popen(
        cwd=cwd,
        stdin=stdin.fileno(),
        stdout=stdout.fileno(),
        stderr=stderr.fileno(),
        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP | DETACHED_PROCESS,
    )
