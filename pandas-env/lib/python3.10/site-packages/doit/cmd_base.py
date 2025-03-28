import inspect
import sys
from collections import deque
from collections import defaultdict
import textwrap

from .globals import Globals
from . import version
from .cmdparse import CmdOption, CmdParse
from .exceptions import InvalidCommand, InvalidDodoFile
from .dependency import CHECKERS, DbmDB, JsonDB, SqliteDB, Dependency, JSONCodec
from .action import CmdAction
from .plugin import PluginDict
from . import loader


def version_tuple(ver_in):
    """convert a version string or tuple into a 3-element tuple with ints
    Any part that is not a number (dev0, a2, b4) will be converted to -1
    """
    result = []
    if isinstance(ver_in, str):
        parts = ver_in.split('.')
    else:
        parts = ver_in
    for rev in parts:
        try:
            result.append(int(rev))
        except ValueError:
            result.append(-1)
    assert len(result) == 3
    return result


def _wrap(content, indent_level):
    """wrap multiple lines keeping the indentation"""
    indent = ' ' * indent_level
    wrap_opt = {
        'initial_indent': indent,
        'subsequent_indent': indent,
    }
    lines = []
    for paragraph in content.splitlines():
        if not paragraph:
            lines.append('')
            continue
        lines.extend(textwrap.wrap(paragraph, **wrap_opt))
    return lines


class Command(object):
    """third-party should subclass this for commands that do no use tasks

    :cvar name: (str) name of sub-cmd to be use from cmdline
    :cvar doc_purpose: (str) single line cmd description
    :cvar doc_usage: (str) describe accepted parameters
    :cvar doc_description: (str) long description/help for cmd
    :cvar cmd_options:
          (list of dict) see cmdparse.CmdOption for dict format
    """

    # if not specified uses the class name
    name = None

    # doc attributes, should be sub-classed
    doc_purpose = ''
    doc_usage = ''
    doc_description = None  # None value will completely omit line from doc

    # sequence of dicts
    cmd_options = tuple()

    # `execute_tasks` indicates whether this command execute task's actions.
    # This is used by the loader to indicate when delayed task creation
    # should be used.
    execute_tasks = False

    def __init__(self, config=None, bin_name='doit', opt_vals=None, **kwargs):
        """configure command

        :param bin_name: str - name of command line program
        :param config: dict

        Set extra configuration values, this vals can come from:
         * directly passed when using the API - through DoitMain.run()
         * from an INI configuration file
        """
        self.bin_name = bin_name
        self.name = self.get_name()
        # config includes all option values and plugins
        self.config = config if config else {}
        self._cmdparser = None
        # option values (i.e. loader options)
        self.opt_vals = opt_vals if opt_vals else {}

        # config_vals contains cmd option values
        self.config_vals = {}
        if 'GLOBAL' in self.config:
            self.config_vals.update(self.config['GLOBAL'])
        if self.name in self.config:
            self.config_vals.update(self.config[self.name])

        # Use post-mortem PDB in case of error loading tasks.
        # Only available for `run` command.
        self.pdb = False


    @classmethod
    def get_name(cls):
        """get command name as used from command line"""
        return cls.name or cls.__name__.lower()

    @property
    def cmdparser(self):
        """get CmdParser instance for this command

        initialize option values:
          - Default are taken from harded option definition
          - Defaults are overwritten from user's cfg (INI) file
        """
        if not self._cmdparser:
            self._cmdparser = CmdParse(self.get_options())
            self._cmdparser.overwrite_defaults(self.config_vals)
        return self._cmdparser


    def get_options(self):
        """@reutrn list of CmdOption
        """
        return [CmdOption(opt) for opt in self.cmd_options]


    def execute(self, opt_values, pos_args):  # pragma: no cover
        """execute command
        :param opt_values: (dict) with cmd_options values
        :param pos_args: (list) of cmd-line positional arguments
        """
        raise NotImplementedError()


    def parse_execute(self, in_args):
        """helper. just parse parameters and execute command

        @args: see method parse
        @returns: result of self.execute
        """
        params, args = self.cmdparser.parse(in_args)
        self.pdb = params.get('pdb', False)
        params.update(self.opt_vals)
        return self.execute(params, args)

    def help(self):
        """return help text"""
        text = []
        text.append("PURPOSE")
        text.extend(_wrap(self.doc_purpose, 4))

        text.append("\nUSAGE")
        usage = "{} {} {}".format(self.bin_name, self.name, self.doc_usage)
        text.extend(_wrap(usage, 4))

        text.append("\nOPTIONS")
        options = defaultdict(list)
        for opt in self.cmdparser.options:
            options[opt.section].append(opt)
        for section, opts in sorted(options.items()):
            section_name = '\n{}'.format(section or self.name)
            text.extend(_wrap(section_name, 4))
            for opt in opts:
                # ignore option that cant be modified on cmd line
                if not (opt.short or opt.long):
                    continue
                text.extend(_wrap(opt.help_param(), 6))
                # TODO It should always display option's default value
                opt_help = opt.help % {'default': opt.default}
                opt_choices = opt.help_choices()
                opt_config = 'config: {}'.format(opt.name)
                if opt.env_var:
                    opt_env = ', environ: {}'.format(opt.env_var)
                else:
                    opt_env = ''
                desc = '{} {} ({}{})'.format(opt_help, opt_choices,
                                             opt_config, opt_env)
                text.extend(_wrap(desc, 12))

                # print bool inverse option
                if opt.inverse:
                    text.extend(_wrap('--{}'.format(opt.inverse), 6))
                    text.extend(_wrap('opposite of --{}'.format(opt.long), 12))

        if self.doc_description is not None:
            text.append("\n\nDESCRIPTION")
            text.extend(_wrap(self.doc_description, 4))
        return "\n".join(text)


######################################################################

# choose internal dependency file.
opt_depfile = {
    'section': 'DB backend',
    'name': 'dep_file',
    'short': '',
    'long': 'db-file',
    'type': str,
    'default': ".doit.db",
    'help': "file used to save successful runs [default: %(default)s]"
}

# dependency file DB backend
opt_backend = {
    'section': 'DB backend',
    'name': 'backend',
    'short': '',
    'long': 'backend',
    'type': str,
    'default': "dbm",
    'help': ("Select dependency file backend. [default: %(default)s]")
}

# dependency file codecs
opt_codec = {
    'section': 'doit core',
    'name': 'codec_cls',
    'short': '',
    'long': '',
    'type': str,
    'default': "json",
    'help': ("Select codec for task's data in database. [default: %(default)s]")
}

opt_check_file_uptodate = {
    'section': 'doit core',
    'name': 'check_file_uptodate',
    'short': '',
    'long': 'check_file_uptodate',
    'type': str,
    'default': 'md5',
    'help': """\
Choose how to check if files have been modified.
Available options [default: %(default)s]:
  'md5': use the md5sum
  'timestamp': use the timestamp
"""
}



#### options related to dodo.py
# select dodo file containing tasks
opt_dodo = {
    'section': 'task loader',
    'name': 'dodoFile',
    'short': 'f',
    'long': 'file',
    'type': str,
    'default': 'dodo.py',
    'env_var': 'DOIT_FILE',
    'help': "load task from dodo FILE [default: %(default)s]"
}

# cwd
opt_cwd = {
    'section': 'task loader',
    'name': 'cwdPath',
    'short': 'd',
    'long': 'dir',
    'type': str,
    'default': None,
    'help': ("set path to be used as cwd directory "
             "(file paths on dodo file are relative to dodo.py location).")
}

# seek dodo file on parent folders
opt_seek_file = {
    'section': 'task loader',
    'name': 'seek_file',
    'short': 'k',
    'long': 'seek-file',
    'type': bool,
    'default': False,
    'env_var': 'DOIT_SEEK_FILE',
    'help': ("seek dodo file on parent folders [default: %(default)s]")
}


class TaskLoader():
    def __init__(self):
        raise NotImplementedError(
            'doit.cmd_base.py:TaskLoader was removed on 0.36.0, use TaskLoader2 instead')

class TaskLoader2():
    """Interface of task loaders with new-style API.

    :cvar cmd_options:
          (list of dict) see cmdparse.CmdOption for dict format

    This API separates the loading of the configuration and the loading
    of the actual tasks, which enables additional elements to be available
    during task creation.
    """
    API = 2
    cmd_options = ()

    def __init__(self):
        # list of command names, used to detect clash of task names and commands
        self.cmd_names = []
        self.config = None   # reference to config object taken from Command
        self.task_opts = None  # dict with task options (no need parsing, API usage)

    def setup(self, opt_values):
        """Delayed initialization.

        To be implemented if the data is needed by derived classes.

        :param opt_values: (dict) with values for cmd_options
        """
        pass

    def load_doit_config(self):
        """Load doit configuration.

        The method must not be called before invocation of ``setup``.

        :return: (dict) Dictionary of doit configuration values.
        """
        raise NotImplementedError()  # pragma: no cover

    def load_tasks(self, cmd, pos_args):
        """Load tasks.

        The method must not be called before invocation of ``load_doit_config``.

        :param cmd: (doit.cmd_base.Command) current command being executed
        :param pos_args: (list str) positional arguments from command line
        :return: (List[Task])
        """
        raise NotImplementedError()  # pragma: no cover


class NamespaceTaskLoader(TaskLoader2):
    """Implementation of a loader of tasks from an abstract namespace.

    A namespace is simply a dictionary to objects like functions and
    objects. See the derived classes for some concrete namespace types.
    """
    def __init__(self):
        super().__init__()
        self.namespace = None

    def load_doit_config(self):
        return loader.load_doit_config(self.namespace)

    def load_tasks(self, cmd, pos_args):
        tasks = loader.load_tasks(
            self.namespace, self.cmd_names, allow_delayed=cmd.execute_tasks,
            args=pos_args, config=self.config, task_opts=self.task_opts)

        # Add task options from config, if present
        if self.config is not None:
            for task in tasks:
                task_stanza = 'task:' + task.name
                if task_stanza in self.config:
                    task.cfg_values = self.config[task_stanza]

        # add values from API run_tasks() usage
        if self.task_opts is not None:
            for task in tasks:
                if self.task_opts and task.name in self.task_opts:
                    task.cfg_values = self.task_opts[task.name]
                    if task.pos_arg and task.pos_arg in task.cfg_values:
                        task.pos_arg_val = task.cfg_values[task.pos_arg]
        return tasks


class ModuleTaskLoader(NamespaceTaskLoader):
    """load tasks from a module/dictionary containing task generators
    Usage: `ModuleTaskLoader(my_module)` or `ModuleTaskLoader(globals())`
    """
    def __init__(self, mod_dict):
        super().__init__()
        if inspect.ismodule(mod_dict):
            self.namespace = dict(inspect.getmembers(mod_dict))
        else:
            self.namespace = mod_dict


class DodoTaskLoader(NamespaceTaskLoader):
    """default task-loader create tasks from a dodo.py file"""
    cmd_options = (opt_dodo, opt_cwd, opt_seek_file)

    def setup(self, opt_values):
        # lazily load namespace from dodo file per config parameters:
        self.namespace = dict(inspect.getmembers(loader.get_module(
            opt_values['dodoFile'],
            opt_values['cwdPath'],
            opt_values['seek_file'],
        )))


def get_loader(config, task_loader=None, cmds=None):
    """get task loader and configure it

    :param config: (dict) the whole config from INI
    :param task_loader: a TaskLoader class
    :param cmds: dict of available commands
    """
    config = config if config else {}
    loader = None
    if task_loader:
        loader = task_loader  # task_loader set from the API
    else:
        global_config = config.get('GLOBAL', {})
        if 'loader' in global_config:
            # a plugin loader
            loader_name = global_config['loader']
            plugins = PluginDict()
            plugins.add_plugins(config, 'LOADER')
            loader = plugins.get_plugin(loader_name)()

    if not loader:
        loader = DodoTaskLoader()  # default loader

    if cmds:
        loader.cmd_names = list(sorted(cmds.keys()))
    loader.config = config
    return loader


################################

class DoitCmdBase(Command):
    """
    subclass must define:
    cmd_options => list of option dictionary (see CmdOption)
    _execute => method, argument names must be option names
    """
    base_options = (opt_depfile, opt_backend, opt_codec,
                    opt_check_file_uptodate)

    def __init__(self, task_loader, cmds=None, **kwargs):
        super(DoitCmdBase, self).__init__(**kwargs)
        self.sel_tasks = None  # selected tasks for command
        self.sel_default_tasks = True  # False if tasks were specified from command line
        self.dep_manager = None
        self.outstream = sys.stdout
        self.loader = task_loader
        self._backends = self.get_backends()


    def get_options(self):
        """from base class - merge base_options, loader_options and cmd_options
        """
        opt_list = (self.base_options + self.loader.cmd_options + self.cmd_options)
        return [CmdOption(opt) for opt in opt_list]


    def _execute(self):  # pragma: no cover
        """to be subclassed - actual command implementation"""
        raise NotImplementedError


    @staticmethod
    def check_minversion(minversion):
        """check if this version of doit satisfy minimum required version
        Minimum version specified by configuration on dodo.
        """
        if minversion:
            if version_tuple(minversion) > version_tuple(version.VERSION):
                msg = ('Please update doit. '
                       'Minimum version required is {required}. '
                       'You are using {actual}. ')
                raise InvalidDodoFile(msg.format(required=minversion,
                                                 actual=version.VERSION))

    def get_checker_cls(self, check_file_uptodate):
        """return checker class to be used by dep_manager"""
        if isinstance(check_file_uptodate, str):
            if check_file_uptodate not in CHECKERS:
                msg = ("No check_file_uptodate named '{}'."
                       " Type '{} help run' to see a list "
                       "of available checkers.").format(
                           check_file_uptodate, self.bin_name)
                raise InvalidCommand(msg)
            return CHECKERS[check_file_uptodate]
        else:
            # user defined class
            return check_file_uptodate

    def get_codec_cls(self, codec):
        """return a class used to encode or decode python-action results"""
        if isinstance(codec, str):
            if codec == 'json':
                return JSONCodec
            else:  # pragma: no cover
                raise NotImplementedError('Implement codec plugin')
        else:
            # user specified class
            return codec


    def get_backends(self):
        """return PluginDict of DB backends, including core and plugins"""
        backend_map = {'dbm': DbmDB, 'json': JsonDB, 'sqlite3': SqliteDB}
        # add plugins
        plugins = PluginDict()
        plugins.add_plugins(self.config, 'BACKEND')
        backend_map.update(plugins.to_dict())

        # set choices, sub-classes might not have this option
        if 'backend' in self.cmdparser:
            choices = {k: getattr(v, 'desc', '') for k, v in backend_map.items()}
            self.cmdparser['backend'].choices = choices

        return backend_map


    def execute(self, params, args):
        """load dodo.py, set attributes and call self._execute

        :param params: instance of cmdparse.DefaultUpdate
        :param args: list of string arguments (containing task names)
        """
        self.loader.setup(params)
        dodo_config = self.loader.load_doit_config()

        # merge config values from dodo.py into params
        params.update_defaults(dodo_config)

        self.check_minversion(params.get('minversion'))

        # set selected tasks for command
        self.sel_default_tasks = len(args) == 0
        self.sel_tasks = args or params.get('default_tasks')

        CmdAction.STRING_FORMAT = params.get('action_string_formatting', 'old')
        if CmdAction.STRING_FORMAT not in ('old', 'both', 'new'):
            raise InvalidDodoFile(
                '`action_string_formatting` must be one of `old`, `both`, `new`')

        # create dep manager
        db_class = self._backends.get(params['backend'])
        checker_cls = self.get_checker_cls(params['check_file_uptodate'])
        codec_cls = self.get_codec_cls(params['codec_cls'])
        # note the command have the responsibility to call dep_manager.close()

        if self.dep_manager is None:
            # dep_manager might have been already set (used on unit-test)
            self.dep_manager = Dependency(
                db_class, params['dep_file'], checker_cls=checker_cls,
                codec_cls=codec_cls)

        # register dependency manager in global registry:
        Globals.dep_manager = self.dep_manager
        # load tasks
        self.task_list = self.loader.load_tasks(cmd=self, pos_args=args)

        # hack to pass parameter into _execute() calls that are not part
        # of command line options
        params['pos_args'] = args
        params['continue_'] = params.get('continue')
        # hack: determine if value came from command line or config
        params['force_verbosity'] = 'verbosity' in params._non_default_keys

        # magic - create dict based on signature of _execute() method.
        # this done so that _execute() have a nice API with name parameters
        # instead of just taking a dict.
        args_name = list(inspect.signature(self._execute).parameters.keys())
        exec_params = dict((n, params[n]) for n in args_name)
        return self._execute(**exec_params)



# helper functions to find list of tasks


def check_tasks_exist(tasks, name_list, skip_wildcard=False):
    """check task exist"""
    if not name_list:
        return
    for task_name in name_list:
        if skip_wildcard and '*' in task_name:
            continue
        if task_name not in tasks:
            msg = "'%s' is not a task."
            raise InvalidCommand(msg % task_name)


# this is used by commands that do not execute tasks (forget, auto...)
def tasks_and_deps_iter(tasks, sel_tasks, yield_duplicates=False):
    """iterator of select_tasks and its dependencies

    @param tasks (dict - Task)
    @param sel_tasks(list - str)
    """
    processed = set()  # str - task name
    to_process = deque(sel_tasks)  # str - task name
    # get initial task
    while to_process:
        task = tasks[to_process.popleft()]
        processed.add(task.name)
        yield task
        # FIXME this does not take calc_dep into account
        for task_dep in task.task_dep + task.setup_tasks:
            if (task_dep not in processed) and (task_dep not in to_process):
                to_process.append(task_dep)
            elif yield_duplicates:
                yield tasks[task_dep]


def subtasks_iter(tasks, task):
    """find all subtasks for a given task
    @param tasks (dict - Task)
    @param task (Task)
    """
    for name in task.task_dep:
        dep = tasks[name]
        if dep.subtask_of == task.name:
            yield dep
