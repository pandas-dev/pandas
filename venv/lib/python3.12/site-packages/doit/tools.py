"""extra goodies to be used in dodo files"""

import os
import time as time_module
import datetime
import json
import hashlib
import operator
import subprocess

from . import exceptions
from .action import CmdAction, PythonAction
from .task import result_dep  # imported for backward compatibility
result_dep  # pyflakes


# action
def create_folder(dir_path):
    """create a folder in the given path if it doesnt exist yet."""
    os.makedirs(dir_path, exist_ok=True)


# title
def title_with_actions(task):
    """return task name task actions"""
    if task.actions:
        title = "\n\t".join([str(action) for action in task.actions])
    # A task that contains no actions at all
    # is used as group task
    else:
        title = "Group: %s" % ", ".join(task.task_dep)
    return "%s => %s" % (task.name, title)



# uptodate
def run_once(task, values):
    """execute task just once
    used when user manually manages a dependency
    """
    def save_executed():
        return {'run-once': True}
    task.value_savers.append(save_executed)
    return values.get('run-once', False)



# uptodate
class config_changed(object):
    """check if passed config was modified
    @var config (str) or (dict)
    @var encoder (json.JSONEncoder) Encoder used to convert non-default values.
    """
    def __init__(self, config, encoder=None):
        self.config = config
        self.config_digest = None
        self.encoder = encoder

    def _calc_digest(self):
        if isinstance(self.config, str):
            return self.config
        elif isinstance(self.config, dict):
            data = json.dumps(self.config, sort_keys=True, cls=self.encoder)
            byte_data = data.encode("utf-8")
            return hashlib.md5(byte_data).hexdigest()
        else:
            msg = ('Invalid type of config_changed parameter got %s,'
                   ' must be string or dict')
            raise Exception(msg % (type(self.config),))

    def configure_task(self, task):
        task.value_savers.append(lambda: {'_config_changed': self.config_digest})

    def __call__(self, task, values):
        """return True if config values are UNCHANGED"""
        self.config_digest = self._calc_digest()
        last_success = values.get('_config_changed')
        if last_success is None:
            return False
        return (last_success == self.config_digest)



# uptodate
class timeout(object):
    """add timeout to task

    @param timeout_limit: (datetime.timedelta, int) in seconds

    if the time elapsed since last time task was executed is bigger than
    the "timeout" time the task is NOT up-to-date
    """

    def __init__(self, timeout_limit):
        if isinstance(timeout_limit, datetime.timedelta):
            self.limit_sec = ((timeout_limit.days * 24 * 3600)
                              + timeout_limit.seconds)
        elif isinstance(timeout_limit, int):
            self.limit_sec = timeout_limit
        else:
            msg = "timeout should be datetime.timedelta or int got %r "
            raise Exception(msg % timeout_limit)

    def __call__(self, task, values):
        def save_now():
            return {'success-time': time_module.time()}
        task.value_savers.append(save_now)
        last_success = values.get('success-time', None)
        if last_success is None:
            return False
        return (time_module.time() - last_success) < self.limit_sec



# uptodate
class check_timestamp_unchanged(object):
    """check if timestamp of a given file/dir is unchanged since last run.

    The C{cmp_op} parameter can be used to customize when timestamps are
    considered unchanged, e.g. you could pass L{operator.ge} to also consider
    e.g. files reverted to an older copy as unchanged; or pass a custom
    function to completely customize what unchanged means.

    If the specified file does not exist, an exception will be raised.  Note
    that if the file C{fn} is a target of another task you should probably add
    C{task_dep} on that task to ensure the file is created before checking it.
    """
    def __init__(self, file_name, time='mtime', cmp_op=operator.eq):
        """initialize the callable

        @param fn: (str) path to file/directory to check
        @param time: (str) which timestamp field to check, can be one of
                     (atime, access, ctime, status, mtime, modify)
        @param cmp_op: (callable) takes two parameters (prev_time, current_time)
                   should return True if the timestamp is considered unchanged

        @raises ValueError: if invalid C{time} value is passed
        """
        if time in ('atime', 'access'):
            self._timeattr = 'st_atime'
        elif time in ('ctime', 'status'):
            self._timeattr = 'st_ctime'
        elif time in ('mtime', 'modify'):
            self._timeattr = 'st_mtime'
        else:
            raise ValueError('time can be one of: atime, access, ctime, '
                             'status, mtime, modify (got: %r)' % time)
        self._file_name = file_name
        self._cmp_op = cmp_op
        self._key = '.'.join([self._file_name, self._timeattr])

    def _get_time(self):
        return getattr(os.stat(self._file_name), self._timeattr)

    def __call__(self, task, values):
        """register action that saves the timestamp and check current timestamp

        @raises OSError: if cannot stat C{self._file_name} file
                         (e.g. doesn't exist)
        """
        def save_now():
            return {self._key: self._get_time()}
        task.value_savers.append(save_now)

        prev_time = values.get(self._key)
        if prev_time is None:  # this is first run
            return False
        current_time = self._get_time()
        return self._cmp_op(prev_time, current_time)


# action class
class LongRunning(CmdAction):
    """Action to handle a Long running shell process,
    usually a server or service.
    Properties:

        * the output is never captured
        * it is always successful (return code is not used)
        * "swallow" KeyboardInterrupt
    """
    def execute(self, out=None, err=None):
        action = self.expand_action()
        process = subprocess.Popen(action, shell=self.shell, **self.pkwargs)
        try:
            process.wait()
        except KeyboardInterrupt:
            # normal way to stop interactive process
            pass

# the name InteractiveAction is deprecated on 0.25
InteractiveAction = LongRunning


class Interactive(CmdAction):
    """Action to handle Interactive shell process:

       * the output is never captured
    """
    def execute(self, out=None, err=None):
        action = self.expand_action()
        process = subprocess.Popen(action, shell=self.shell, **self.pkwargs)
        process.wait()
        if process.returncode != 0:
            return exceptions.TaskFailed(
                "Interactive command failed: '%s' returned %s" %
                (action, process.returncode))



# action class
class PythonInteractiveAction(PythonAction):
    """Action to handle Interactive python:

       * the output is never captured
       * it is successful unless a exception is raised
    """
    def execute(self, out=None, err=None):
        kwargs = self._prepare_kwargs()
        try:
            returned_value = self.py_callable(*self.args, **kwargs)
        except Exception as exception:
            return exceptions.TaskError("PythonAction Error", exception)
        if isinstance(returned_value, str):
            self.result = returned_value
        elif isinstance(returned_value, dict):
            self.values = returned_value
            self.result = returned_value


# debug helper
def set_trace():  # pragma: no cover
    """start debugger, make sure stdout shows pdb output.
    output is not restored.
    """
    import pdb
    import sys
    debugger = pdb.Pdb(stdin=sys.__stdin__, stdout=sys.__stdout__)
    debugger.set_trace(sys._getframe().f_back)  # pylint: disable=W0212



def load_ipython_extension(ip=None):  # pragma: no cover
    """
    Defines a ``%doit`` magic function[1] that discovers and execute tasks
    from IPython's interactive variables (global namespace).

    It will fail if not invoked from within an interactive IPython shell.

    .. Tip::
        To permanently add this magic-function to your IPython, create a new
        script inside your startup-profile
        (``~/.ipython/profile_default/startup/doit_magic.ipy``) with the
        following content:

            %load_ext doit
            %reload_ext doit
            %doit list

    [1] http://ipython.org/ipython-doc/dev/interactive/tutorial.html#magic-functions
    """
    from IPython.core.getipython import get_ipython
    from IPython.core.magic import register_line_magic

    from doit.cmd_base import ModuleTaskLoader
    from doit.doit_cmd import DoitMain

    # Only (re)load_ext provides the ip context.
    ip = ip or get_ipython()

    @register_line_magic
    def doit(line):
        """
        Run *doit* with `task_creators` from all interactive variables
        (IPython's global namespace).

        Examples:

            >>> %doit --help          ## Show help for options and arguments.

            >>> def task_foo():
                    return {'actions': ['echo hi IPython'],
                            'verbosity': 2}

            >>> %doit list            ## List any tasks discovered.
            foo

            >>> %doit                 ## Run any tasks.
            .  foo
            hi IPython

        """
        # Override db-files location inside ipython-profile dir,
        # which is certainly writable.
        prof_dir = ip.profile_dir.location
        opt_vals = {'dep_file': os.path.join(prof_dir, 'db', '.doit.db')}
        commander = DoitMain(ModuleTaskLoader(ip.user_module),
                             extra_config={'GLOBAL': opt_vals})
        commander.BIN_NAME = 'doit'
        commander.run(line.split())

# also expose another way of registering ipython extension
register_doit_as_IPython_magic = load_ipython_extension
