# Licensed under a 3-clause BSD style license - see LICENSE.rst

class ProfilerGui:
    """
    A base class to define a Profiler GUI that is available through
    the ``asv profile`` command.
    """
    name = None
    description = None

    @classmethod
    def is_available(cls):
        """
        Return `True` is the given GUI tool appears to be installed
        and available.
        """
        raise NotImplementedError()

    @classmethod
    def open_profiler_gui(cls, profiler_file):
        """
        Open the profiler GUI to display the results in the given
        profiler file.
        """
        raise NotImplementedError()
