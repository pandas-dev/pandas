"""doit - Automation Tool

The MIT License

Copyright (c) 2008-present Eduardo Naufel Schettino

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from doit.version import VERSION

__version__ = VERSION


from doit import loader
from doit.loader import create_after, task_params
from doit.doit_cmd import get_var
from doit.api import run
from doit.tools import load_ipython_extension
from doit.globals import Globals


__all__ = ['get_var', 'run', 'create_after', 'task_params', 'Globals']

def get_initial_workdir():
    """working-directory from where the doit command was invoked on shell"""
    return loader.initial_workdir

assert load_ipython_extension  # silence pyflakes
