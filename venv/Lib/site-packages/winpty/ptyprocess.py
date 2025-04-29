# -*- coding: utf-8 -*-

# Standard library imports
import codecs
import os
import shlex
import signal
import socket
import subprocess
import threading
import time
from shutil import which

# Local imports
from .winpty import PTY


class PtyProcess(object):
    """This class represents a process running in a pseudoterminal.

    The main constructor is the :meth:`spawn` classmethod.
    """

    def __init__(self, pty):
        assert isinstance(pty, PTY)
        self.pty = pty
        self.pid = pty.pid
        # self.fd = pty.fd

        self.read_blocking = bool(int(os.environ.get('PYWINPTY_BLOCK', 1)))
        self.closed = False
        self.flag_eof = False

        # Used by terminate() to give kernel time to update process status.
        # Time in seconds.
        self.delayafterterminate = 0.1
        # Used by close() to give kernel time to update process status.
        # Time in seconds.
        self.delayafterclose = 0.1

        # Set up our file reader sockets.
        self._server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._server.bind(("127.0.0.1", 0))
        address = self._server.getsockname()
        self._server.listen(1)

        # Read from the pty in a thread.
        self._thread = threading.Thread(target=_read_in_thread,
            args=(address, self.pty, self.read_blocking))
        self._thread.daemon = True
        self._thread.start()

        self.fileobj, _ = self._server.accept()
        self.fd = self.fileobj.fileno()

    @classmethod
    def spawn(cls, argv, cwd=None, env=None, dimensions=(24, 80),
              backend=None):
        """Start the given command in a child process in a pseudo terminal.

        This does all the setting up the pty, and returns an instance of
        PtyProcess.

        Dimensions of the psuedoterminal used for the subprocess can be
        specified as a tuple (rows, cols), or the default (24, 80) will be
        used.
        """
        if isinstance(argv, str):
            argv = shlex.split(argv, posix=False)

        if not isinstance(argv, (list, tuple)):
            raise TypeError("Expected a list or tuple for argv, got %r" % argv)

        # Shallow copy of argv so we can modify it
        argv = argv[:]
        command = argv[0]
        env = env or os.environ

        path = env.get('PATH', os.defpath)
        command_with_path = which(command, path=path)
        if command_with_path is None:
            raise FileNotFoundError(
                'The command was not found or was not ' +
                'executable: %s.' % command
            )
        command = command_with_path
        argv[0] = command
        cmdline = ' ' + subprocess.list2cmdline(argv[1:])
        cwd = cwd or os.getcwd()

        backend = backend or os.environ.get('PYWINPTY_BACKEND', None)
        backend = int(backend) if backend is not None else backend

        proc = PTY(dimensions[1], dimensions[0],
                   backend=backend)

        # Create the environemnt string.
        envStrs = []
        for (key, value) in env.items():
            envStrs.append('%s=%s' % (key, value))
        env = '\0'.join(envStrs) + '\0'

        # command = bytes(command, encoding)
        # cwd = bytes(cwd, encoding)
        # cmdline = bytes(cmdline, encoding)
        # env = bytes(env, encoding)

        if len(argv) == 1:
            proc.spawn(command, cwd=cwd, env=env)
        else:
            proc.spawn(command, cwd=cwd, env=env, cmdline=cmdline)

        inst = cls(proc)
        inst._winsize = dimensions

        # Set some informational attributes
        inst.argv = argv
        if env is not None:
            inst.env = env
        if cwd is not None:
            inst.launch_dir = cwd

        return inst

    @property
    def exitstatus(self):
        """The exit status of the process.
        """
        return self.pty.get_exitstatus()

    def fileno(self):
        """This returns the file descriptor of the pty for the child.
        """
        return self.fd

    def close(self, force=False):
        """This closes the connection with the child application. Note that
        calling close() more than once is valid. This emulates standard Python
        behavior with files. Set force to True if you want to make sure that
        the child is terminated (SIGKILL is sent if the child ignores
        SIGINT)."""
        if not self.closed:
            self.fileobj.close()
            self._server.close()
            # Give kernel time to update process status.
            time.sleep(self.delayafterclose)
            if self.isalive():
                if not self.terminate(force):
                    raise IOError('Could not terminate the child.')
            self.fd = -1
            self.closed = True
            # del self.pty

    def __del__(self):
        """This makes sure that no system resources are left open. Python only
        garbage collects Python objects. OS file descriptors are not Python
        objects, so they must be handled explicitly. If the child file
        descriptor was opened outside of this class (passed to the constructor)
        then this does not close it.
        """
        # It is possible for __del__ methods to execute during the
        # teardown of the Python VM itself. Thus self.close() may
        # trigger an exception because os.close may be None.
        try:
            self.close()
        except Exception:
            pass

    def flush(self):
        """This does nothing. It is here to support the interface for a
        File-like object. """
        pass

    def isatty(self):
        """This returns True if the file descriptor is open and connected to a
        tty(-like) device, else False."""
        return self.isalive()

    def read(self, size=1024):
        """Read and return at most ``size`` characters from the pty.

        Can block if there is nothing to read. Raises :exc:`EOFError` if the
        terminal was closed.
        """
        # try:
        #     data = self.pty.read(size, blocking=self.read_blocking)
        # except Exception as e:
        #     if "EOF" in str(e):
        #         raise EOFError(e) from e
        # return data
        data = self.fileobj.recv(size)
        if not data:
            self.flag_eof = True
            raise EOFError('Pty is closed')

        if data == b'0011Ignore':
            data = ''

        err = True
        while err and data:
            try:
                data.decode('utf-8')
                err = False
            except UnicodeDecodeError:
                data += self.fileobj.recv(1)
        return data.decode('utf-8')

    def readline(self):
        """Read one line from the pseudoterminal as bytes.

        Can block if there is nothing to read. Raises :exc:`EOFError` if the
        terminal was closed.
        """
        buf = []
        while 1:
            try:
                ch = self.read(1)
            except EOFError:
                return ''.join(buf)
            buf.append(ch)
            if ch == '\n':
                return ''.join(buf)

    def write(self, s):
        """Write the string ``s`` to the pseudoterminal.

        Returns the number of bytes written.
        """
        if not self.pty.isalive():
            raise EOFError('Pty is closed')

        nbytes = self.pty.write(s)
        return nbytes

    def terminate(self, force=False):
        """This forces a child process to terminate."""
        if not self.isalive():
            return True
        self.kill(signal.SIGINT)
        time.sleep(self.delayafterterminate)
        if not self.isalive():
            return True
        if force:
            self.kill(signal.SIGTERM)
            time.sleep(self.delayafterterminate)
            if not self.isalive():
                return True
            else:
                return False

    def wait(self):
        """This waits until the child exits. This is a blocking call. This will
        not read any data from the child.
        """
        while self.isalive():
            time.sleep(0.1)
        return self.exitstatus

    def isalive(self):
        """This tests if the child process is running or not. This is
        non-blocking. If the child was terminated then this will read the
        exitstatus or signalstatus of the child. This returns True if the child
        process appears to be running or False if not.
        """
        alive = self.pty.isalive()
        self.closed = not alive
        return alive

    def kill(self, sig=None):
        """Kill the process with the given signal.
        """
        os.kill(self.pid, sig)

    def sendcontrol(self, char):
        '''Helper method that wraps send() with mnemonic access for sending control
        character to the child (such as Ctrl-C or Ctrl-D).  For example, to send
        Ctrl-G (ASCII 7, bell, '\a')::
            child.sendcontrol('g')
        See also, sendintr() and sendeof().
        '''
        char = char.lower()
        a = ord(char)
        if 97 <= a <= 122:
            a = a - ord('a') + 1
            byte = bytes([a]).decode("ascii")
            return self.pty.write(byte), byte
        d = {'@': 0, '`': 0,
            '[': 27, '{': 27,
            '\\': 28, '|': 28,
            ']': 29, '}': 29,
            '^': 30, '~': 30,
            '_': 31,
            '?': 127}
        if char not in d:
            return 0, ''

        byte = bytes([d[char]]).decode("ascii")
        return self.pty.write(byte), byte

    def sendeof(self):
        """This sends an EOF to the child. This sends a character which causes
        the pending parent output buffer to be sent to the waiting child
        program without waiting for end-of-line. If it is the first character
        of the line, the read() in the user program returns 0, which signifies
        end-of-file. This means to work as expected a sendeof() has to be
        called at the beginning of a line. This method does not send a newline.
        It is the responsibility of the caller to ensure the eof is sent at the
        beginning of a line."""
        # Send control character 4 (Ctrl-D)
        self.pty.write('\x04'), '\x04'

    def sendintr(self):
        """This sends a SIGINT to the child. It does not require
        the SIGINT to be the first character on a line. """
        # Send control character 3 (Ctrl-C)
        self.pty.write('\x03'), '\x03'

    def eof(self):
        """This returns True if the EOF exception was ever raised.
        """
        return self.flag_eof

    def getwinsize(self):
        """Return the window size of the pseudoterminal as a tuple (rows, cols).
        """
        return self._winsize

    def setwinsize(self, rows, cols):
        """Set the terminal window size of the child tty.
        """
        self._winsize = (rows, cols)
        self.pty.set_size(cols, rows)


def _read_in_thread(address, pty, blocking):
    """Read data from the pty in a thread.
    """
    client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client.connect(address)

    call = 0

    while 1:
        try:
            data = pty.read(4096, blocking=blocking) or b'0011Ignore'
            try:
                client.send(bytes(data, 'utf-8'))
            except socket.error:
                break

            # Handle end of file.
            if pty.iseof():
                try:
                    client.send(b'')
                except socket.error:
                    pass
                finally:
                    break

            call += 1
        except Exception as e:
            break
        time.sleep(1e-3)

    client.close()
