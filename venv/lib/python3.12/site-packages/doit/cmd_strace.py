import sys
import os
import re

from .exceptions import InvalidCommand
from .action import CmdAction
from .task import Task
from .cmd_run import Run


# filter to display only files from cwd
opt_show_all = {
    'name': 'show_all',
    'short': 'a',
    'long': 'all',
    'type': bool,
    'default': False,
    'help': "display all files (not only from within CWD path)",
}

opt_keep_trace = {
    'name': 'keep_trace',
    'long': 'keep',
    'type': bool,
    'default': False,
    'help': "save strace command output into strace.txt",
}


class Strace(Run):
    doc_purpose = "use strace to list file_deps and targets"
    doc_usage = "TASK"
    doc_description = """
The output is a list of files prefixed with 'R' for open in read mode
or 'W' for open in write mode.
The files are listed in chronological order.

This is a debugging feature with many limitations.
  * can strace only one task at a time
  * can only strace CmdAction
  * the process being traced itself might have some kind of cache,
    that means it might not write a target file if it exist
  * does not handle chdir

So this is NOT 100% reliable, use with care!
"""

    cmd_options = (opt_show_all, opt_keep_trace)

    TRACE_CMD = "strace -f -e trace=file %s 2>>%s "
    TRACE_OUT = 'strace.txt'

    def execute(self, params, args):
        """remove existing output file if any and do sanity checking"""
        if os.path.exists(self.TRACE_OUT):  # pragma: no cover
            os.unlink(self.TRACE_OUT)
        if len(args) != 1:
            msg = ('strace failed, must select *one* task to strace.'
                   '\nCheck `{} help strace`.'.format(self.bin_name))
            raise InvalidCommand(msg)
        result = Run.execute(self, params, args)
        if (not params['keep_trace']) and os.path.exists(self.TRACE_OUT):
            os.unlink(self.TRACE_OUT)
        return result

    def _execute(self, show_all):
        """1) wrap the original action with strace and save output in file
           2) add a second task that will generate the report from temp file
        """
        # find task to trace and wrap it
        selected = self.sel_tasks[0]
        for task in self.task_list:
            if task.name == selected:
                self.wrap_strace(task)
                break

        # add task to print report
        report_strace = Task(
            'strace_report',
            actions=[(find_deps, [self.outstream, self.TRACE_OUT, show_all])],
            verbosity=2,
            task_dep=[selected],
            uptodate=[False],
        )
        self.task_list.append(report_strace)
        self.sel_tasks.append(report_strace.name)

        # clear strace file
        return Run._execute(self, sys.stdout)


    @classmethod
    def wrap_strace(cls, task):
        """wrap task actions into strace command"""
        wrapped_actions = []
        for action in task.actions:
            if isinstance(action, CmdAction):
                cmd = cls.TRACE_CMD % (action._action, cls.TRACE_OUT)
                wrapped = CmdAction(cmd, task, save_out=action.save_out)
                wrapped_actions.append(wrapped)
            else:
                wrapped_actions.append(action)
        task._action_instances = wrapped_actions
        # task should be always executed
        task._extend_uptodate([False])


def find_deps(outstream, strace_out, show_all):
    """read file witn strace output, return dict with deps, targets"""
    # 7978  open("/etc/ld.so.cache", O_RDONLY|O_CLOEXEC) = 3
    # get "mode" file was open, until ')' is closed
    # ignore rest of line
    # .*\(                 # ignore text until '('
    # [^"]*"               # ignore text until '"'
    # (?P<file>[^"]*)"     # get "file" name inside "
    # , (\[.*\])*          # ignore elements if inside [] - used by execve
    # (?P<mode>[^)]*)\)    # get mode opening file
    #  = ].*               # check syscall was successful""",
    regex = re.compile(
        r'.*\([^"]*"(?P<file>[^"]*)",'
        + r' (\[.*\])*(?P<mode>[^)]*)\) = [^-].*')

    read = set()
    write = set()
    cwd = os.getcwd()
    if not os.path.exists(strace_out):
        return
    with open(strace_out) as text:
        for line in text:
            # ignore non file operation
            match = regex.match(line)
            if not match:
                continue
            rel_name = match.group('file')
            name = os.path.abspath(rel_name)

            # ignore files out of cwd
            if not show_all:
                if not name.startswith(cwd):
                    continue

            if 'WR' in match.group('mode'):
                if name not in write:
                    write.add(name)
                    outstream.write("W %s\n" % name)
            else:
                if name not in read:
                    read.add(name)
                    outstream.write("R %s\n" % name)
