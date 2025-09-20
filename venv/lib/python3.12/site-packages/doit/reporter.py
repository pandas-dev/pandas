"""Reports doit execution status/results"""

import sys
import time
import datetime
import json
from io import StringIO

from .exceptions import BaseFail


class ConsoleReporter(object):
    """Default reporter. print results on console/terminal (stdout/stderr)

    @ivar failure_verbosity: (int) include captured stdout/stderr on failure
                             report even if already shown.
    """
    # short description, used by the help system
    desc = 'console output'

    def __init__(self, outstream, options):
        # save non-successful result information (include task errors)
        self.failures = []
        self.runtime_errors = []
        self.failure_verbosity = options.get('failure_verbosity', 0)
        self.outstream = outstream

    def write(self, text):
        self.outstream.write(text)

    def initialize(self, tasks, selected_tasks):
        """called just after tasks have been loaded before execution starts"""
        pass


    def get_status(self, task):
        """called when task is selected (check if up-to-date)"""
        pass

    def execute_task(self, task):
        """called when execution starts"""
        # ignore tasks that do not define actions
        # ignore private/hidden tasks (tasks that start with an underscore)
        if task.actions and (task.name[0] != '_'):
            self.write('.  %s\n' % task.title())

    def add_failure(self, task, fail: BaseFail):
        """called when execution finishes with a failure"""
        result = {'task': task, 'exception': fail}
        if fail.report:
            self.failures.append(result)
            self._write_failure(result)

    def add_success(self, task):
        """called when execution finishes successfully"""
        pass

    def skip_uptodate(self, task):
        """skipped up-to-date task"""
        if task.name[0] != '_':
            self.write("-- %s\n" % task.title())

    def skip_ignore(self, task):
        """skipped ignored task"""
        self.write("!! %s\n" % task.title())

    def cleanup_error(self, exception):
        """error during cleanup"""
        sys.stderr.write(exception.get_msg())

    def runtime_error(self, msg):
        """error from doit (not from a task execution)"""
        # saved so they are displayed after task failures messages
        self.runtime_errors.append(msg)

    def teardown_task(self, task):
        """called when starts the execution of teardown action"""
        pass


    def _write_failure(self, result, write_exception=True):
        msg = '%s - taskid:%s\n' % (result['exception'].get_name(),
                                    result['task'].name)
        self.write(msg)
        if write_exception:
            self.write(result['exception'].get_msg())
            self.write("\n")

    def complete_run(self):
        """called when finished running all tasks"""
        # if test fails print output from failed task
        for result in self.failures:
            task = result['task']
            # makes no sense to print output if task was not executed
            if not task.executed:
                continue
            show_err = task.verbosity < 1 or self.failure_verbosity > 0
            show_out = task.verbosity < 2 or self.failure_verbosity == 2
            if show_err or show_out:
                self.write("#" * 40 + "\n")
            if show_err:
                self._write_failure(result,
                                    write_exception=self.failure_verbosity)
                err = "".join([a.err for a in task.actions if a.err])
                self.write("{} <stderr>:\n{}\n".format(task.name, err))
            if show_out:
                out = "".join([a.out for a in task.actions if a.out])
                self.write("{} <stdout>:\n{}\n".format(task.name, out))

        if self.runtime_errors:
            self.write("#" * 40 + "\n")
            self.write("Execution aborted.\n")
            self.write("\n".join(self.runtime_errors))
            self.write("\n")


class ExecutedOnlyReporter(ConsoleReporter):
    """No output for skipped (up-to-date) and group tasks

    Produces zero output unless a task is executed
    """
    desc = 'console, no output for skipped (up-to-date) and group tasks'

    def skip_uptodate(self, task):
        """skipped up-to-date task"""
        pass

    def skip_ignore(self, task):
        """skipped ignored task"""
        pass



class ZeroReporter(ConsoleReporter):
    """Report only internal errors from doit"""
    desc = 'report only internal errors from doit'

    def _just_pass(self, *args):
        """over-write base to do nothing"""
        pass

    get_status = execute_task = add_failure = add_success \
        = skip_uptodate = skip_ignore = teardown_task = complete_run \
        = _just_pass

    def runtime_error(self, msg):
        sys.stderr.write(msg)


class ErrorOnlyReporter(ZeroReporter):
    desc = """Report only errors internal or TaskError and TaskFailure."""

    def add_failure(self, task, fail_info: BaseFail):
        if not fail_info.report:
            return
        exception_name = fail_info.get_name()
        self.write(f'taskid:{task.name} - {exception_name}\n')
        self.write(fail_info.get_msg())
        self.write("\n")


class TaskResult(object):
    """result object used by JsonReporter"""
    # FIXME what about returned value from python-actions ?
    def __init__(self, task):
        self.task = task
        self.result = None  # fail, success, up-to-date, ignore
        self.out = None  # stdout from task
        self.err = None  # stderr from task
        self.error = None  # error from doit (exception traceback)
        self.started = None  # datetime when task execution started
        self.elapsed = None  # time (in secs) taken to execute task
        self._started_on = None  # timestamp
        self._finished_on = None  # timestamp

    def start(self):
        """called when task starts its execution"""
        self._started_on = time.time()

    def set_result(self, result, error=None):
        """called when task finishes its execution"""
        self._finished_on = time.time()
        self.result = result
        line_sep = "\n<------------------------------------------------>\n"
        self.out = line_sep.join([a.out for a in self.task.actions if a.out])
        self.err = line_sep.join([a.err for a in self.task.actions if a.err])
        self.error = error

    def to_dict(self):
        """convert result data to dictionary"""
        if self._started_on is not None:
            started = datetime.datetime.utcfromtimestamp(self._started_on)
            self.started = str(started.strftime('%Y-%m-%d %H:%M:%S.%f'))
            self.elapsed = self._finished_on - self._started_on
        return {'name': self.task.name,
                'result': self.result,
                'out': self.out,
                'err': self.err,
                'error': self.error,
                'started': self.started,
                'elapsed': self.elapsed}


class JsonReporter(object):
    """output results in JSON format

    - out (str)
    - err (str)
    - tasks (list - dict):
         - name (str)
         - result (str)
         - out (str)
         - err (str)
         - error (str)
         - started (str)
         - elapsed (float)
    """

    desc = 'output in JSON format'

    def __init__(self, outstream, options=None):  # pylint: disable=W0613
        # options parameter is not used
        # json result is sent to stdout when doit finishes running
        self.t_results = {}
        # when using json reporter output can not contain any other output
        # than the json data. so anything that is sent to stdout/err needs to
        # be captured.
        self._old_out = sys.stdout
        sys.stdout = StringIO()
        self._old_err = sys.stderr
        sys.stderr = StringIO()
        self.outstream = outstream
        # runtime and cleanup errors
        self.errors = []

    def get_status(self, task):
        """called when task is selected (check if up-to-date)"""
        self.t_results[task.name] = TaskResult(task)

    def execute_task(self, task):
        """called when execution starts"""
        self.t_results[task.name].start()

    def add_failure(self, task, exception):
        """called when execution finishes with a failure"""
        self.t_results[task.name].set_result('fail', exception.get_msg())

    def add_success(self, task):
        """called when execution finishes successfully"""
        self.t_results[task.name].set_result('success')

    def skip_uptodate(self, task):
        """skipped up-to-date task"""
        self.t_results[task.name].set_result('up-to-date')

    def skip_ignore(self, task):
        """skipped ignored task"""
        self.t_results[task.name].set_result('ignore')

    def cleanup_error(self, exception):
        """error during cleanup"""
        self.errors.append(exception.get_msg())

    def runtime_error(self, msg):
        """error from doit (not from a task execution)"""
        self.errors.append(msg)

    def teardown_task(self, task):
        """called when starts the execution of teardown action"""
        pass

    def complete_run(self):
        """called when finished running all tasks"""
        # restore stdout
        log_out = sys.stdout.getvalue()
        sys.stdout = self._old_out
        log_err = sys.stderr.getvalue()
        sys.stderr = self._old_err

        # add errors together with stderr output
        if self.errors:
            log_err += "\n".join(self.errors)

        task_result_list = [
            tr.to_dict() for tr in self.t_results.values()]
        json_data = {'tasks': task_result_list,
                     'out': log_out,
                     'err': log_err}
        # indent not available on simplejson 1.3 (debian etch)
        # json.dump(json_data, sys.stdout, indent=4)
        json.dump(json_data, self.outstream)
