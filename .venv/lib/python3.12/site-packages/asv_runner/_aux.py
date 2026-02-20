import contextlib
import importlib
import os
import sys
import tempfile

from .benchmarks._maxrss import set_cpu_affinity


class SpecificImporter:
    """
    Module importer that only allows loading a given module from the
    given path.

    #### Notes
    Using this enables importing the asv benchmark suite without
    adding its parent directory to sys.path. The parent directory can
    in principle contain anything, including some version of the
    project module (common situation if asv.conf.json is on project
    repository top level).
    """

    def __init__(self, name, root):
        """
        Initialize a new instance of `SpecificImporter`.

        #### Parameters
        **name** (`str`)
        : The name of the module to load.

        **root** (`str`)
        : The path to the directory containing the module.
        """
        self._name = name
        self._root = root

    def find_spec(self, fullname, path, target):
        """
        Find the module specification for the given module.

        #### Parameters
        **fullname** (`str`)
        : The fully qualified name of the module.

        **path** (list or None)
        : The path for module search, or None if unavailable.

        **target** (object)
        : The target object to import.

        #### Returns
        **spec** (`ModuleSpec` or None)
        : The module specification if the module is found, or None otherwise.

        #### Notes
        This method is called by the import system to find the module
        specification for the requested module.  If the requested module matches
        the name of the SpecificImporter instance, it returns the module
        specification using the `importlib.machinery.PathFinder`.
        """
        if fullname == self._name:
            if path is not None:
                raise ValueError()
            finder = importlib.machinery.PathFinder()
            return finder.find_spec(fullname, [self._root], target)
        return None


def update_sys_path(root):
    """
    Update sys.meta_path to include the SpecificImporter.

    ##### Parameters
     `root` (`str`): The path to the root directory.

    ##### Notes
    This function inserts the SpecificImporter into the `sys.meta_path` at the
    beginning, allowing the module to be imported using the SpecificImporter
    when it is encountered during the import process.
    """
    sys.meta_path.insert(
        0, SpecificImporter(os.path.basename(root), os.path.dirname(root))
    )


@contextlib.contextmanager
def posix_redirect_output(filename=None, permanent=True):
    """
    Redirect stdout/stderr to a file, using posix `dup2`.

    #### Parameters
    **filename** (`str` or None, optional)
    : The name of the file to redirect the output to. If None, a temporary
      file will be created.

    **permanent** (`bool`, optional)
    : Indicates whether the redirection is permanent or temporary. If False,
      the original stdout/stderr will be restored after the context is exited.

    #### Yields
    **filename** (`str`)
    : The name of the file where the output is redirected.

    #### Notes
    The function redirects the `stdout` and `stderr` streams to a file using
    the posix `dup2` function. It is typically used within a `with` statement to
    encapsulate the code block where the redirection is desired.

    If `filename` is not provided, a temporary file will be created and used for
    redirection.

    If `permanent` is `True`, the redirection will persist after the context is
    exited. If `False`, the original `stdout`/`stderr` will be restored.
    """
    sys.stdout.flush()
    sys.stderr.flush()

    stdout_fd = sys.stdout.fileno()
    stderr_fd = sys.stderr.fileno()

    if not permanent:
        stdout_fd_copy = os.dup(stdout_fd)
        stderr_fd_copy = os.dup(stderr_fd)

    if filename is None:
        out_fd, filename = tempfile.mkstemp()
    else:
        out_fd = os.open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)

    try:
        # Redirect stdout and stderr to file
        os.dup2(out_fd, stdout_fd)
        os.dup2(out_fd, stderr_fd)

        yield filename
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.close(out_fd)

        if not permanent:
            os.dup2(stdout_fd_copy, stdout_fd)
            os.dup2(stderr_fd_copy, stderr_fd)
            os.close(stdout_fd_copy)
            os.close(stderr_fd_copy)


def recvall(sock, size):
    """
    Receive data of given size from a socket connection.

    #### Parameters
    **sock** (socket object)
    : The socket connection to receive data from.

    **size** (`int`)
    : The size of the data to receive, in bytes.

    #### Returns
    **data** (`bytes`)
    : The received data.

    #### Raises
    **RuntimeError**
    : If the data received from the socket is less than the specified size.

    #### Notes
    The function receives data from a socket connection in multiple chunks until
    the specified size is reached. It ensures that all the required data is received
    before returning.

    If the received data size is less than the specified size, a `RuntimeError`
    is raised indicating the failure to receive the complete data.
    """
    data = b""
    while len(data) < size:
        s = sock.recv(size - len(data))
        data += s
        if not s:
            raise RuntimeError(
                "did not receive data from socket " f"(size {size}, got only {data!r})"
            )
    return data


def set_cpu_affinity_from_params(extra_params):
    """
    Set CPU affinity based on the provided parameters.

    #### Parameters
    **extra_params** (`dict` or `None`)
    : Additional parameters containing CPU affinity information.

    #### Notes
    This function attempts to set the CPU affinity for the current process
    based on the provided parameters. It uses the `set_cpu_affinity` function
    internally to perform the actual affinity setting.

    If the `extra_params` dictionary contains a key "cpu_affinity" with a
    valid affinity list, the CPU affinity will be set accordingly.

    #### Raises
    **BaseException**
    : If setting the CPU affinity fails, an exception is raised and an error
      message is printed.

    #### Example
    ```{code-block} python
    extra_params = {"cpu_affinity": [0, 1]}
    set_cpu_affinity_from_params(extra_params)
    ```
    """
    affinity_list = extra_params.get("cpu_affinity", None)
    if affinity_list is not None:
        try:
            set_cpu_affinity(affinity_list)
        except BaseException as exc:
            print(f"asv: setting cpu affinity {affinity_list !r} failed: {exc !r}")
