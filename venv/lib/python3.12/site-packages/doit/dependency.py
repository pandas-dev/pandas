"""Manage (save/check) task dependency-on-files data."""

import os
import hashlib
import subprocess
import inspect
import json
from collections import defaultdict
from dbm import dumb
import dbm as ddbm

# uncomment imports below to run tests on all dbm backends...
# import dumbdbm as ddbm
# import dbm as ddbm
# import gdbm as ddbm

# note: to check which DBM backend is being used (in py2):
#       >>> anydbm._defaultmod



class DatabaseException(Exception):
    """Exception class for whatever backend exception"""
    pass


def get_md5(input_data):
    """return md5 from string or unicode"""
    byte_data = input_data.encode("utf-8")
    return hashlib.md5(byte_data).hexdigest()


def get_file_md5(path):
    """Calculate the md5 sum from file content.

    @param path: (string) file path
    @return: (string) md5
    """
    with open(path, 'rb') as file_data:
        md5 = hashlib.md5()
        block_size = 128 * md5.block_size
        while True:
            data = file_data.read(block_size)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


class JSONCodec():
    """default implmentation for codec used to save individual task's data"""
    def __init__(self):
        self.encoder = json.JSONEncoder()
        self.decoder = json.JSONDecoder()

    def encode(self, data):
        return self.encoder.encode(data)

    def decode(self, data):
        return self.decoder.decode(data)



class JsonDB(object):
    """Backend using a single text file with JSON content"""

    def __init__(self, name, codec):
        """Open/create a DB file"""
        self.name = name
        self.codec = codec
        if not os.path.exists(self.name):
            self._db = {}
        else:
            self._db = self._load()

    def _load(self):
        """load db content from file"""
        db_file = open(self.name, 'r')
        try:
            try:
                return self.codec.decode(db_file.read())
            except ValueError as error:
                # file contains corrupted json data
                fname = os.path.abspath(self.name)
                msg = (f"{error.args[0]}\nInvalid JSON data in {fname}\n"
                       "To fix this problem, you can just remove the "
                       "corrupted file, a new one will be generated.\n")
                error.args = (msg,)
                raise DatabaseException(msg)
        finally:
            db_file.close()

    def dump(self):
        """save DB content in file"""
        try:
            db_file = open(self.name, 'w')
            db_file.write(self.codec.encode(self._db))
        finally:
            db_file.close()

    def set(self, task_id, dependency, value):
        """Store value in the DB."""
        if task_id not in self._db:
            self._db[task_id] = {}
        self._db[task_id][dependency] = value


    def get(self, task_id, dependency):
        """Get value stored in the DB.

        @return: (string) or (None) if entry not found
        """
        if task_id in self._db:
            return self._db[task_id].get(dependency, None)


    def in_(self, task_id):
        """@return bool if task_id is in DB"""
        return task_id in self._db


    def remove(self, task_id):
        """remove saved dependencies from DB for taskId"""
        if task_id in self._db:
            del self._db[task_id]

    def remove_all(self):
        """remove saved dependencies from DB for all tasks"""
        self._db = {}



class DbmDB(object):
    """Backend using a DBM file with individual values encoded in JSON

    On initialization all items are read from DBM file and loaded on ``_dbm``.
    During execution whenever an item is read (``get`` method) the `json` value
    is cached on ``_db``.
    If an item is modified ``_db`` is update and the `id` is added
    to the `dirty` set. Only on ``dump`` all dirty items values are encoded
    in json into ``_dbm`` and the DBM file is saved.

    :ivar str name: file name/path
    :ivar dbm _dbm: items with json encoded values
    :ivar dict _db: items with python-dict as value
    :ivar set dirty: id of modified tasks
    """
    DBM_CONTENT_ERROR_MSG = 'db type could not be determined'

    def __init__(self, name, codec):
        """Open/create a DB file"""
        self.name = name
        self.codec = codec
        try:
            self._dbm = ddbm.open(self.name, 'c')
        except ddbm.error as exception:
            message = str(exception)
            if message == self.DBM_CONTENT_ERROR_MSG:
                # When a corrupted/old format database is found
                # suggest the user to just remove the file
                new_message = (
                    'Dependencies file in %(filename)s seems to use '
                    'an old format or is corrupted.\n'
                    'To fix the issue you can just remove the database file(s) '
                    'and a new one will be generated.'
                    % {'filename': repr(self.name)})
                raise DatabaseException(new_message)
            else:
                # Re-raise any other exceptions
                raise DatabaseException(message)

        self._db = {}
        self.dirty = set()

    def dump(self):
        """save/close DBM file"""
        for task_id in self.dirty:
            self._dbm[task_id] = self.codec.encode(self._db[task_id])
        self._dbm.close()


    def set(self, task_id, dependency, value):
        """Store value in the DB."""
        if task_id not in self._db:
            self._db[task_id] = {}
        self._db[task_id][dependency] = value
        self.dirty.add(task_id)


    def _in_dbm(self, key):
        """
        should be just::
          return key in self._dbm

         for get()/set() key is convert to bytes but not for 'in'
        """
        return key.encode('utf-8') in self._dbm


    def get(self, task_id, dependency):
        """Get value stored in the DB.

        :return: string or None if entry not found
        """
        # optimization, just try to get it without checking it exists
        if task_id in self._db:
            return self._db[task_id].get(dependency, None)
        else:
            try:
                task_data = self._dbm[task_id]
            except KeyError:
                return
            self._db[task_id] = self.codec.decode(task_data.decode('utf-8'))
            return self._db[task_id].get(dependency, None)


    def in_(self, task_id):
        """@return bool if task_id is in DB"""
        return self._in_dbm(task_id) or task_id in self.dirty


    def remove(self, task_id):
        """remove saved dependencies from DB for taskId"""
        if task_id in self._db:
            del self._db[task_id]
        if self._in_dbm(task_id):
            del self._dbm[task_id]
        if task_id in self.dirty:
            self.dirty.remove(task_id)


    def remove_all(self):
        """remove saved dependencies from DB for all tasks"""
        self._db = {}
        # dumb dbm always opens file in update mode
        if isinstance(self._dbm, dumb._Database):  # pragma: no cover
            self._dbm._index = {}
            self._dbm.close()
        # gdbm can not be running on 2 instances on same thread
        # see https://bitbucket.org/schettino72/doit/issue/16/
        del self._dbm
        self._dbm = ddbm.open(self.name, 'n')
        self.dirty = set()



class SqliteDB(object):
    """ sqlite3 json backend """

    def __init__(self, name, codec):
        self.name = name
        self.codec = codec
        self._conn = self._sqlite3(self.name)
        self._cache = {}
        self._dirty = set()

    def _sqlite3(self, name):
        """Open/create a sqlite3 DB file"""

        # Import sqlite here so it's only imported when required
        import sqlite3
        def dict_factory(cursor, row):
            """convert row to dict"""
            data = {}
            for idx, col in enumerate(cursor.description):
                data[col[0]] = row[idx]
            return data
        def converter(data):
            return self.codec.decode(data.decode('utf-8'))

        sqlite3.register_adapter(list, self.codec.encode)
        sqlite3.register_adapter(dict, self.codec.encode)
        sqlite3.register_converter("json", converter)
        conn = sqlite3.connect(
            name,
            detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
            isolation_level='DEFERRED')
        conn.row_factory = dict_factory
        sqlscript = """
            create table if not exists doit (
                task_id text not null primary key,
                task_data json
            );"""
        try:
            conn.execute(sqlscript)
        except sqlite3.DatabaseError as exception:
            new_message = (
                'Dependencies file in %(filename)s seems to use '
                'an bad format or is corrupted.\n'
                'To fix the issue you can just remove the database file(s) '
                'and a new one will be generated.'
                'Original error: %(msg)s'
                % {'filename': repr(name), 'msg': str(exception)})
            raise DatabaseException(new_message)
        return conn

    def get(self, task_id, dependency):
        """Get value stored in the DB.

        @return: (string) or (None) if entry not found
        """
        if task_id in self._cache:
            return self._cache[task_id].get(dependency, None)
        else:
            data = self._cache[task_id] = self._get_task_data(task_id)
            return data.get(dependency, None)

    def _get_task_data(self, task_id):
        data = self._conn.execute('select task_data from doit where task_id=?',
                                  (task_id,)).fetchone()
        return data['task_data'] if data else {}

    def set(self, task_id, dependency, value):
        """Store value in the DB."""
        if task_id not in self._cache:
            self._cache[task_id] = {}
        self._cache[task_id][dependency] = value
        self._dirty.add(task_id)


    def in_(self, task_id):
        if task_id in self._cache:
            return True
        if self._conn.execute('select task_id from doit where task_id=?',
                              (task_id,)).fetchone():
            return True
        return False

    def dump(self):
        """save/close sqlite3 DB file"""
        for task_id in self._dirty:
            self._conn.execute('insert or replace into doit values (?,?)',
                               (task_id, self.codec.encode(self._cache[task_id])))
        self._conn.commit()
        self._conn.close()
        self._dirty = set()

    def remove(self, task_id):
        """remove saved dependencies from DB for taskId"""
        if task_id in self._cache:
            del self._cache[task_id]
        if task_id in self._dirty:
            self._dirty.remove(task_id)
        self._conn.execute('delete from doit where task_id=?', (task_id,))

    def remove_all(self):
        """remove saved dependencies from DB for all task"""
        self._conn.execute('delete from doit')
        self._cache = {}
        self._dirty = set()


class FileChangedChecker(object):
    """Base checker for dependencies, must be inherited."""

    CheckerError = os.error

    def exists(self, file_path):
        return os.path.exists(file_path)

    def info(self, file_path):
        return os.stat(file_path)

    def check_modified(self, file_path, file_stat, state):
        """Check if file in file_path is modified from previous "state".

        @param file_path (string): file path
        @param file_stat: result of os.stat() of file_path
        @param state: state that was previously saved with ``get_state()``
        @returns (bool): True if dep is modified

        """
        raise NotImplementedError()

    def get_state(self, dep, current_state):
        """Compute the state of a task after it has been successfully executed.

        @param dep (str): path of the dependency file.
        @param current_state (tuple): the current state, saved from a previous
            execution of the task (None if the task was never run).
        @returns: the new state. Return None if the state is unchanged.

        The parameter `current_state` is passed to allow speed optimization,
        see MD5Checker.get_state().
        """
        raise NotImplementedError()


class MD5Checker(FileChangedChecker):
    """MD5 checker, uses the md5sum.

    This is the default checker used by doit.

    As an optimization the check uses (timestamp, file-size, md5).
    If the timestamp is the same it considers that the file has the same
    content. If file size is different its content certainly is modified.
    Finally the md5 is used for a different timestamp with the same size.
    """

    def check_modified(self, file_path, file_stat, state):
        """Check if file in file_path is modified from previous "state".
        """
        timestamp, size, file_md5 = state

        # 1 - if timestamp is not modified file is the same
        if file_stat.st_mtime == timestamp:
            return False

        # 2 - if size is different file is modified
        if file_stat.st_size != size:
            return True

        # 3 - check md5
        return file_md5 != get_file_md5(file_path)


    def get_state(self, dep, current_state):
        timestamp = os.path.getmtime(dep)
        # time optimization. if dep is already saved with current
        # timestamp skip calculating md5
        if current_state and current_state[0] == timestamp:
            return
        size = os.path.getsize(dep)
        md5 = get_file_md5(dep)
        return timestamp, size, md5


class TimestampChecker(FileChangedChecker):
    """Checker that use only the timestamp."""

    def check_modified(self, file_path, file_stat, state):
        return file_stat.st_mtime != state

    def get_state(self, dep, current_state):
        """@returns float: mtime for file `dep`"""
        return os.path.getmtime(dep)


# name of checkers class available
CHECKERS = {'md5': MD5Checker,
            'timestamp': TimestampChecker}


class DependencyStatus(object):
    """Result object for Dependency.get_status.

    @ivar status: (str) one of "run", "up-to-date" or "error"
    """

    def __init__(self, get_log):
        self.get_log = get_log
        self.status = 'up-to-date'
        # save reason task is not up-to-date
        self.reasons = defaultdict(list)
        self.error_reason = None

    def add_reason(self, reason, arg, status='run'):
        """sets state and append reason for not being up-to-date
        :return boolean: processing should be interrupted
        """
        self.status = status
        if self.get_log:
            self.reasons[reason].append(arg)
        return not self.get_log

    def set_reason(self, reason, arg):
        """sets state and reason for not being up-to-date
        :return boolean: processing should be interrupted
        """
        self.status = 'run'
        if self.get_log:
            self.reasons[reason] = arg
        return not self.get_log

    def get_error_message(self):
        '''return str with error message'''
        return self.error_reason



class Dependency(object):
    """Manage tasks dependencies.

    Each dependency is saved in "db". There are several "db" backends.
    It uses a Key-Value format where the key is task-name
    and value is a dictionary.
    Each task has a dictionary where keys are `dependency`'s (absolute file path),
    and the value is the dependency signature.
    Apart from dependencies other values are also saved on the task dictionary:

    * ``_values_:`` task's values
    * ``result:`` task result

    And also some internal doit attributes:

    * ``ignore:``
    * ``deps:``
    * ``checker:``

    Those can be accessed with generic DB ``get()``, see below...

    :ivar string name: filepath of the DB file
    :ivar bool _closed: DB was flushed to file
    """
    def __init__(self, db_class, backend_name, checker_cls=MD5Checker,
                 codec_cls=JSONCodec):
        self._closed = False
        self.checker = checker_cls()
        self.db_class = db_class
        self.backend = db_class(backend_name, codec=codec_cls())
        self._set = self.backend.set
        self._get = self.backend.get
        self.remove = self.backend.remove
        self.remove_all = self.backend.remove_all
        self._in = self.backend.in_
        self.name = self.backend.name

    def close(self):
        """Write DB in file"""
        if not self._closed:
            self.backend.dump()
            self._closed = True


    ####### task specific

    def save_success(self, task, result_hash=None):
        """save info after a task is successfully executed

        :param str result_hash: explicitly set result_hash
        """
        # save task values
        self._set(task.name, "_values_:", task.values)

        # save task result md5
        if result_hash is not None:
            self._set(task.name, "result:", result_hash)
        elif task.result:
            if isinstance(task.result, dict):
                self._set(task.name, "result:", task.result)
            else:
                self._set(task.name, "result:", get_md5(task.result))

        # file-dep
        self._set(task.name, 'checker:', self.checker.__class__.__name__)
        for dep in task.file_dep:
            state = self.checker.get_state(dep, self._get(task.name, dep))
            if state is not None:
                self._set(task.name, dep, state)

        # save list of file_deps
        self._set(task.name, 'deps:', tuple(task.file_dep))

    def get_values(self, task_name):
        """get all saved values from a task

        :return dict:
        """
        values = self._get(task_name, '_values_:')
        return values or {}

    def get_value(self, task_id, key_name):
        """get saved value from task

        :param str task_id:
        :param str key_name: key result dict of the value
        """
        if not self._in(task_id):
            # FIXME do not use generic exception
            raise Exception("taskid '%s' has no computed value!" % task_id)
        values = self.get_values(task_id)
        if key_name not in values:
            msg = "Invalid arg name. Task '%s' has no value for '%s'."
            raise Exception(msg % (task_id, key_name))
        return values[key_name]

    def get_result(self, task_name):
        """get the result saved from a task

        :return (dict or md5sum):
        """
        return self._get(task_name, 'result:')

    def remove_success(self, task):
        """remove saved info from task"""
        self.remove(task.name)

    def ignore(self, task):
        """mark task to be ignored"""
        self._set(task.name, 'ignore:', '1')

    def status_is_ignore(self, task):
        """check if task is marked to be ignored"""
        return self._get(task.name, "ignore:")

    def get_status(self, task, tasks_dict, get_log=False):
        """Check if task is up to date. set task.dep_changed

        If the checker class changed since the previous run, the task is
        deleted, to be sure that its state is not re-used.

        @param task: (Task)
        @param tasks_dict: (dict: Task) passed to objects used on uptodate
        @param get_log: (bool) if True, adds all reasons to the return
                               object why this file will be rebuild.
        @return: (DependencyStatus) a status object with possible status
                                    values up-to-date, run or error

        task.dep_changed (list-strings): file-dependencies that are not
        up-to-date if task not up-to-date because of a target, returned value
        will contain all file-dependencies regardless they are up-to-date
        or not.
        """
        result = DependencyStatus(get_log)
        task.dep_changed = []

        # check uptodate bool/callables
        uptodate_result_list = []
        for utd, utd_args, utd_kwargs in task.uptodate:
            # if parameter is a callable
            if hasattr(utd, '__call__'):
                # FIXME control verbosity, check error messages
                # 1) setup object with global info all tasks
                if isinstance(utd, UptodateCalculator):
                    utd.setup(self, tasks_dict)
                # 2) add magic positional args for `task` and `values`
                # if present.
                spec_args = list(inspect.signature(utd).parameters.keys())
                magic_args = []
                for i, name in enumerate(spec_args):
                    if i == 0 and name == 'task':
                        magic_args.append(task)
                    elif i == 1 and name == 'values':
                        magic_args.append(self.get_values(task.name))
                args = magic_args + utd_args
                # 3) call it and get result
                uptodate_result = utd(*args, **utd_kwargs)
            elif isinstance(utd, str):
                uptodate_result = subprocess.call(
                    utd, shell=True,
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL) == 0
            # parameter is a value
            else:
                uptodate_result = utd

            # None means uptodate was not really calculated and should be
            # just ignored
            if uptodate_result is None:
                continue
            uptodate_result_list.append(uptodate_result)
            if not uptodate_result:
                result.add_reason('uptodate_false', (utd, utd_args, utd_kwargs))

        # any uptodate check is false
        if not get_log and result.status == 'run':
            return result

        # no dependencies means it is never up to date.
        if not (task.file_dep or uptodate_result_list):
            if result.set_reason('has_no_dependencies', True):
                return result


        # if target file is not there, task is not up to date
        for targ in task.targets:
            if not self.checker.exists(targ):
                task.dep_changed = list(task.file_dep)
                if result.add_reason('missing_target', targ):
                    return result

        # check for modified file_dep checker
        previous = self._get(task.name, 'checker:')
        checker_name = self.checker.__class__.__name__
        if previous and previous != checker_name:
            task.dep_changed = list(task.file_dep)
            # remove all saved values otherwise they might be re-used by
            # some optimization on MD5Checker.get_state()
            self.remove(task.name)
            if result.set_reason('checker_changed', (previous, checker_name)):
                return result

        # check for modified file_dep
        previous = self._get(task.name, 'deps:')
        previous_set = set(previous) if previous else None
        if previous_set and previous_set != task.file_dep:
            if get_log:
                added_files = sorted(list(task.file_dep - previous_set))
                removed_files = sorted(list(previous_set - task.file_dep))
                result.set_reason('added_file_dep', added_files)
                result.set_reason('removed_file_dep', removed_files)
            result.status = 'run'

        # list of file_dep that changed
        check_modified = self.checker.check_modified
        changed = []
        for dep in task.file_dep:
            state = self._get(task.name, dep)
            try:
                file_stat = self.checker.info(dep)
            except self.checker.CheckerError:
                error_msg = "Dependent file '{}' does not exist.".format(dep)
                result.error_reason = error_msg.format(dep)
                if result.add_reason('missing_file_dep', dep, 'error'):
                    return result
            else:
                if state is None or check_modified(dep, file_stat, state):
                    changed.append(dep)
        task.dep_changed = changed

        if len(changed) > 0:
            result.set_reason('changed_file_dep', changed)

        return result



#############

class UptodateCalculator(object):
    """Base class for 'uptodate' that need access to all tasks
    """
    def __init__(self):
        self.get_val = None  # Dependency._get
        self.tasks_dict = None  # dict with all tasks

    def setup(self, dep_manager, tasks_dict):
        """@param"""
        self.get_val = dep_manager._get
        self.tasks_dict = tasks_dict
