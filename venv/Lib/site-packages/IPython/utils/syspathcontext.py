import sys
import warnings


class prepended_to_syspath:
    """A context for prepending a directory to sys.path for a second."""

    def __init__(self, dir):
        self.dir = dir

    def __enter__(self):
        if self.dir not in sys.path:
            sys.path.insert(0, self.dir)
            self.added = True
        else:
            self.added = False

    def __exit__(self, type, value, traceback):
        if self.added:
            try:
                sys.path.remove(self.dir)
            except ValueError:
                pass
        # Returning False causes any exceptions to be re-raised.
        return False
