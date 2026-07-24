# Licensed under a 3-clause BSD style license - see LICENSE.rst

from .. import asv_profiling, util


class RunSnakeRunGui(asv_profiling.ProfilerGui):
    name = 'runsnake'
    description = "RunSnakeRun http://www.vrplumber.com/programming/runsnakerun/"

    @classmethod
    def is_available(cls):
        return util.has_command('runsnake')

    @classmethod
    def open_profiler_gui(cls, profiler_file):
        command = util.which('runsnake')

        return util.check_call([command, profiler_file], timeout=None)
