"""Loads dodo file (a python module) and convert them to 'tasks' """

import os
import sys
import copy
import inspect
import importlib
from collections import OrderedDict

from .exceptions import InvalidTask, InvalidCommand, InvalidDodoFile
from .task import DelayedLoader, Task, dict_to_task
from .cmdparse import TaskParse, CmdOption


# Directory path from where doit was executed.
# Set by loader, to be used on dodo.py by users.
initial_workdir = None

# TASK_STRING: (string) prefix used to identify python function
# that are task generators in a dodo file.
TASK_STRING = "task_"

def flat_generator(gen, gen_doc=''):
    """return only values from generators
    if any generator yields another generator it is recursively called
    """
    for item in gen:
        if inspect.isgenerator(item):
            item_doc = item.gi_code.co_consts[0]
            for value, value_doc in flat_generator(item, item_doc):
                yield value, value_doc
        else:
            yield item, gen_doc



def get_module(dodo_file, cwd=None, seek_parent=False):
    """
    Find python module defining tasks, it is called "dodo" file.

    @param dodo_file(str): path to file containing the tasks
    @param cwd(str): path to be used cwd, if None use path from dodo_file
    @param seek_parent(bool): search for dodo_file in parent paths if not found
    @return (module) dodo module
    """
    global initial_workdir
    initial_workdir = os.getcwd()
    def exist_or_raise(path):
        """raise exception if file on given path doesnt exist"""
        if not os.path.exists(path):
            msg = (f"Could not find dodo file '{path}'.\n"
                   "Please use '-f' to specify file name.\n")
            raise InvalidDodoFile(msg)

    # get absolute path name
    if os.path.isabs(dodo_file):
        dodo_path = dodo_file
        exist_or_raise(dodo_path)
    else:
        if not seek_parent:
            dodo_path = os.path.abspath(dodo_file)
            exist_or_raise(dodo_path)
        else:
            # try to find file in any folder above
            current_dir = initial_workdir
            dodo_path = os.path.join(current_dir, dodo_file)
            file_name = os.path.basename(dodo_path)
            parent = os.path.dirname(dodo_path)
            while not os.path.exists(dodo_path):
                new_parent = os.path.dirname(parent)
                if new_parent == parent:  # reached root path
                    exist_or_raise(dodo_file)
                parent = new_parent
                dodo_path = os.path.join(parent, file_name)

    ## load module dodo file and set environment
    base_path, file_name = os.path.split(dodo_path)
    # make sure dodo path is on sys.path so we can import it
    sys.path.insert(0, base_path)

    if cwd is None:
        # by default cwd is same as dodo.py base path
        full_cwd = base_path
    else:
        # insert specified cwd into sys.path
        full_cwd = os.path.abspath(cwd)
        if not os.path.isdir(full_cwd):
            msg = "Specified 'dir' path must be a directory.\nGot '%s'(%s)."
            raise InvalidCommand(msg % (cwd, full_cwd))
        sys.path.insert(0, full_cwd)

    # file specified on dodo file are relative to cwd
    os.chdir(full_cwd)

    # get module containing the tasks
    return importlib.import_module(os.path.splitext(file_name)[0])



def create_after(executed=None, target_regex=None, creates=None):
    """Annotate a task-creator function with delayed loader info"""
    def decorated(func):

        func.doit_create_after = DelayedLoader(
            func,
            executed=executed,
            target_regex=target_regex,
            creates=creates
        )
        return func
    return decorated


def task_params(param_def=None):
    """Annotate a task-creator function with definition of parameters
    to get arguments from cmd line
    """
    if param_def is None or type(param_def) != list:
        raise ValueError('task_params must be called with a valid parameter definition.')
    def decorated(func):
        func._task_creator_params = param_def
        return func
    return decorated



def load_tasks(namespace, command_names=(), allow_delayed=False, args=(),
               config=None, task_opts=None):
    """Find task-creators and create tasks

    @param namespace: (dict) containing the task creators, it might
                        contain other stuff
    @param command_names: (list - str) blacklist for task names
    @param allow_delayed: (bool) if True ignore doit_crate_after['executed']
    @param args: (list - str) command line arguments (task names and option arguments)
    @param config: (dict) configuration taken from TOML and INI files

    `load_all == False` is used by the runner to delay the creation of
    tasks until a dependent task is executed. This is only used by the `run`
    command, other commands should always load all tasks since it wont execute
    any task.

    @return task_list (list) of Tasks in the order they were defined on the file
    """
    funcs = _get_task_creators(namespace, command_names)
    # sort by the order functions were defined (line number)
    # TODO: this ordering doesnt make sense when generators come
    # from different modules
    funcs.sort(key=lambda obj: obj[2])

    task_list = []

    def _append_params(tasks, param_def):
        """Apply parameters defined for the task generator to the tasks
        defined by the generator.
        """
        for task in tasks:
            if task.subtask_of is None:  # only parent tasks
                # task.params can not be used with creator_params
                if task.params:
                    msg = (f"Task '{task.name}'. `params` attribute can not be used"
                           " in conjuction with `@task_params`")
                    raise InvalidTask(msg)
                task.creator_params = param_def


    def _process_gen(ref, creator_kwargs):
        """process a task creator, generating tasks"""
        gen_tasks = generate_tasks(name, ref(**creator_kwargs), ref.__doc__)
        if hasattr(ref, '_task_creator_params'):
            _append_params(gen_tasks, ref._task_creator_params)
        task_list.extend(gen_tasks)

    def _add_delayed(tname, ref, original_delayed, kwargs):
        # Make sure create_after can be used on class methods.
        # delayed.creator is initially set by the decorator,
        # so always an unbound function.
        # Here we re-assign with the reference taken on doit load phase
        # because it is bounded method.
        this_delayed = copy.copy(delayed)
        this_delayed.creator = ref
        d_task = Task(tname, None, loader=this_delayed, doc=delayed.creator.__doc__)

        if hasattr(ref, '_task_creator_params'):
            this_delayed.kwargs = kwargs
            d_task.creator_params = getattr(ref, '_task_creator_params', None)
        task_list.append(d_task)


    # Map arg_name to its position.
    # Save only args that do not start with `-` (potentially task names)
    arg_pos = {}
    for index, term in enumerate(args):
        if term[0] != '-':
            arg_pos[term] = index

    for name, ref, _ in funcs:
        delayed = getattr(ref, 'doit_create_after', None)

        # Parse command line arguments for task generator parameters
        creator_params = getattr(ref, '_task_creator_params', None)
        if creator_params is not None:
            parser = TaskParse([CmdOption(opt) for opt in creator_params])
            # Add task options from config, if present
            if config:
                task_stanza = 'task:' + name
                if task_stanza in config:
                    parser.overwrite_defaults(config[task_stanza])

            # option params passed through API
            if task_opts and name in task_opts:
                creator_kwargs = task_opts[name]
            # if relevant command line defaults are available parse those
            elif name in arg_pos:
                creator_kwargs, _ = parser.parse(args[arg_pos[name] + 1:])
            else:
                creator_kwargs, _ = parser.parse('')
        else:
            creator_kwargs = {}

        if not delayed:  # not a delayed task, just run creator
            _process_gen(ref, creator_kwargs)
        elif delayed.creates:  # delayed with explicit task basename
            for tname in delayed.creates:
                _add_delayed(tname, ref, delayed, creator_kwargs)
        elif allow_delayed:  # delayed no explicit name, cmd run
            _add_delayed(name, ref, delayed, creator_kwargs)
        else:  # delayed no explicit name, cmd list (run creator)
            _process_gen(ref, creator_kwargs)

    return task_list


def _get_task_creators(namespace, command_names):
    """get functions defined in the `namespace` and select the task-creators

    A task-creator is a function that:
       - name starts with the string TASK_STRING
       - has the attribute `create_doit_tasks`

    @return (list - func) task-creators
    """
    funcs = []
    prefix_len = len(TASK_STRING)
    # get all functions that are task-creators
    for name, ref in namespace.items():

        # Do not solicit tasks from the @task_params decorator.
        if ref is task_params:
            continue

        # function is a task creator because of its name
        if name.startswith(TASK_STRING) and (
                inspect.isfunction(ref) or inspect.ismethod(ref)):
            # remove TASK_STRING prefix from name
            task_name = name[prefix_len:]

        # object is a task creator because it contains the special method
        elif hasattr(ref, 'create_doit_tasks'):
            ref = ref.create_doit_tasks
            # create_doit_tasks might have a basename to overwrite task name.
            task_name = getattr(ref, 'basename', name)

        # ignore functions that are not a task creator
        else:  # pragma: no cover
            # coverage can't get "else: continue"
            continue

        # tasks can't have the same name of a commands
        if task_name in command_names:
            msg = (f"Task can't be called '{task_name}' because this is a command name."
                   " Please choose another name.")
            raise InvalidDodoFile(msg)
        # get line number where function is defined
        line = inspect.getsourcelines(ref)[1]
        # add to list task generator functions
        funcs.append((task_name, ref, line))

    return funcs


def load_doit_config(dodo_module):
    """
    @param dodo_module (dict) dict with module members
    """
    doit_config = dodo_module.get('DOIT_CONFIG', {})
    if not isinstance(doit_config, dict):
        msg = ("DOIT_CONFIG  must be a dict. got:'%s'%s")
        raise InvalidDodoFile(msg % (repr(doit_config), type(doit_config)))
    return doit_config


def _generate_task_from_return(func_name, task_dict, gen_doc):
    """generate a single task from a dict returned by a task generator"""
    if 'name' in task_dict:
        raise InvalidTask("Task '%s'. Only subtasks use field name." %
                          func_name)

    task_dict['name'] = task_dict.pop('basename', func_name)

    # Use task generator docstring
    # if no doc present in task dict
    if 'doc' not in task_dict:
        task_dict['doc'] = gen_doc

    return dict_to_task(task_dict)


def _generate_task_from_yield(tasks, func_name, task_dict, gen_doc):
    """generate a single task from a dict yielded by task generator

    @param tasks: dictionary with created tasks
    @return None: the created task is added to 'tasks' dict
    """
    # check valid input
    if not isinstance(task_dict, dict):
        raise InvalidTask("Task '%s' must yield dictionaries" %
                          func_name)

    msg_dup = "Task generation '%s' has duplicated definition of '%s'"
    basename = task_dict.pop('basename', None)
    # if has 'name' this is a sub-task
    if 'name' in task_dict:
        basename = basename or func_name
        # if subname is None attributes from group task
        if task_dict['name'] is None:
            task_dict['name'] = basename
            task_dict['actions'] = None
            group_task = dict_to_task(task_dict)
            group_task.has_subtask = True
            tasks[basename] = group_task
            return

        # name is '<task>.<subtask>'
        full_name = f"{basename}:{task_dict['name']}"
        if full_name in tasks:
            raise InvalidTask(msg_dup % (func_name, full_name))
        task_dict['name'] = full_name
        sub_task = dict_to_task(task_dict)
        sub_task.subtask_of = basename

        # get/create task group
        group_task = tasks.get(basename)
        if group_task:
            if not group_task.has_subtask:
                raise InvalidTask(msg_dup % (func_name, basename))
        else:
            group_task = Task(basename, None, doc=gen_doc, has_subtask=True)
            tasks[basename] = group_task
        group_task.task_dep.append(sub_task.name)
        tasks[sub_task.name] = sub_task
    # NOT a sub-task
    else:
        if not basename:
            raise InvalidTask(
                "Task '%s' must contain field 'name' or 'basename'. %s" %
                (func_name, task_dict))
        if basename in tasks:
            raise InvalidTask(msg_dup % (func_name, basename))
        task_dict['name'] = basename
        # Use task generator docstring if no doc present in task dict
        if 'doc' not in task_dict:
            task_dict['doc'] = gen_doc
        tasks[basename] = dict_to_task(task_dict)


def generate_tasks(func_name, gen_result, gen_doc=None):
    """Create tasks from a task generator result.

    @param func_name: (string) name of taskgen function
    @param gen_result: value returned by a task generator function
                       it can be a dict or generator (generating dicts)
    @param gen_doc: (string/None) docstring from the task generator function
    @param param_def: (dict) additional task parameter definitions
                      passed down from generator
    @return: (list - Task)
    """
    # a task instance, just return it without any processing
    if isinstance(gen_result, Task):
        return (gen_result,)

    # task described as a dictionary
    if isinstance(gen_result, dict):
        return [_generate_task_from_return(func_name, gen_result, gen_doc)]

    # a generator
    if inspect.isgenerator(gen_result):
        tasks = OrderedDict()  # task_name: task
        # the generator return subtasks as dictionaries
        for task_dict, x_doc in flat_generator(gen_result, gen_doc):
            if isinstance(task_dict, Task):
                tasks[task_dict.name] = task_dict
            else:
                _generate_task_from_yield(tasks, func_name, task_dict, x_doc)

        if tasks:
            return list(tasks.values())
        else:
            # special case task_generator did not generate any task
            # create an empty group task
            return [Task(func_name, None, doc=gen_doc, has_subtask=True)]

    if gen_result is None:
        return ()

    raise InvalidTask(
        "Task '%s'. Must return a dictionary or generator. Got %s" %
        (func_name, type(gen_result)))
