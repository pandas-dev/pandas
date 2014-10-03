from __future__ import division, absolute_import, print_function

from setuptools.command.egg_info import egg_info as _egg_info

class egg_info(_egg_info):
    def run(self):
        # We need to ensure that build_src has been executed in order to give
        # setuptools' egg_info command real filenames instead of functions which
        # generate files.
        self.run_command("build_src")
        _egg_info.run(self)
