"""doit CLI (command line interface)"""

import os
import sys
import traceback
from collections import defaultdict
from configparser import ConfigParser
import importlib

from .version import VERSION
from .plugin import PluginDict
from .exceptions import InvalidDodoFile, InvalidCommand, InvalidTask
from .cmdparse import CmdOption, CmdParseError, CmdParse
from .cmd_base import get_loader
from .cmd_help import Help
from .cmd_run import Run
from .cmd_clean import Clean
from .cmd_list import List
from .cmd_info import Info
from .cmd_forget import Forget
from .cmd_ignore import Ignore
from .cmd_dumpdb import DumpDB
from .cmd_strace import Strace
from .cmd_completion import TabCompletion
from .cmd_resetdep import ResetDep


# used to save variable values passed from command line
_CMDLINE_VARS = None

def reset_vars():
    global _CMDLINE_VARS
    _CMDLINE_VARS = {}

def get_var(name, default=None):
    # Ignore if not initialized.
    # This is a work-around for Windows multi-processing
    # See https://github.com/pydoit/doit/issues/164
    if _CMDLINE_VARS is None:
        return None
    return _CMDLINE_VARS.get(name, default)

def set_var(name, value):
    _CMDLINE_VARS[name] = value



class DoitConfig():
    """Parse and store values taken from INI and TOML configuration files"""
    # support TOML python libs
    _TOML_LIBS = ['tomllib', 'tomli', 'tomlkit']
    PLUGIN_TYPES = ['command', 'loader', 'backend', 'reporter']

    def __init__(self):
        self._toml = None
        self.config = defaultdict(dict)

    def loads(self, config_filenames):
        for config_filename in config_filenames:
            if str(config_filename).lower().endswith('.toml'):
                prefix = 'tool.doit' if config_filename == 'pyproject.toml' else ''
                toml_config = self.load_config_toml(config_filename, prefix)
                for section in toml_config:
                    self.config[section].update(toml_config[section].items())
            else:
                # INI config format
                ini_config = self.load_config_ini(config_filename)
                for section in ini_config.sections():
                    self.config[section].update(ini_config[section].items())

    @property
    def toml(self):
        """get available toml lib, if any"""
        if self._toml is None:
            for toml_lib in self._TOML_LIBS:
                try:
                    self._toml = importlib.import_module(toml_lib)
                    break
                except ImportError:
                    pass
        return self._toml

    def as_dict(self):
        """return values in legacy dict format"""
        # FIXME: convert INI to TOML format instead of converting TOML to INI
        return self.config


    @staticmethod
    def load_config_ini(filenames):
        """read config from INI files

        :param files: str or list of str.
                      Like ConfigParser.read() param filenames
        """
        cfg_parser = ConfigParser(allow_no_value=True, delimiters=('=',))
        cfg_parser.optionxform = str  # preserve case of option names
        cfg_parser.read(filenames)
        return cfg_parser


    def load_config_toml(self, filename, prefix):
        """read config from a TOML file, adapt to ConfigParser structure

        :param filename: str

        returns an empty dictionary if tool.doit is not present in the parsed data
        """

        toml_config = {}
        raw = None

        if os.path.exists(filename):
            if not self.toml:
                sys.stderr.write(
                    f'''WARNING: File "{filename}" might contain doit configuration,'''
                    '''but a TOML parser is not available.\n'''
                    f'''\tPlease install one of: {', '.join(self._TOML_LIBS)}.\n'''
                )

            else:
                with open(filename, encoding='utf-8') as fp:
                    text = fp.read()

                raw = self.toml.loads(text)

                if raw and isinstance(raw, dict):
                    if prefix:
                        for part in prefix.split('.'):
                            raw = raw.get(part, {})
                    doit_toml = raw

                    # hoist /tool/doit/plugins/AAA to /AAA
                    for plugin_type, plugins in doit_toml.pop('plugins', {}).items():
                        assert plugin_type in self.PLUGIN_TYPES
                        toml_config[plugin_type.upper()] = plugins

                    # hoist /tool/doit/commands/bbb to /bbb
                    for command, command_config in doit_toml.pop('commands', {}).items():
                        toml_config[command] = command_config

                    # hoist /tool/doit/tasks/ccc to /task:ccc
                    for task, task_config in doit_toml.pop('tasks', {}).items():
                        toml_config['task:{}'.format(task)] = task_config

                    # hoist /tool/doit/ddd to /GLOBAL/ddd
                    for global_name, global_value in doit_toml.items():
                        toml_config.setdefault('GLOBAL', {})[global_name] = global_value

        return toml_config


class DoitMain(object):
    # core doit commands
    BIN_NAME = sys.argv[0].split('/')[-1]
    DOIT_CMDS = (Help, Run, List, Info, Clean, Forget, Ignore, DumpDB,
                 Strace, TabCompletion, ResetDep)

    def __init__(self, task_loader=None,
                 config_filenames=('pyproject.toml', 'doit.cfg'),
                 extra_config=None):
        """
        :param extra_config: dict of extra argument values (by argument name)
                             This is parameter is only used by explicit API call.
        """
        self.task_loader = task_loader

        # backward compability: convert single filename to list
        if isinstance(config_filenames, str):
            config_filenames = [config_filenames]

        # ignore config files do that not exist
        config_filenames = [fn for fn in config_filenames if os.path.exists(fn)]

        self.config = defaultdict(dict)
        if extra_config:
            for section, items in extra_config.items():
                self.config[section].update(items)

        # combine config option from INI/TOML files and API
        config_in = DoitConfig()
        config_in.loads(config_filenames)
        for section, vals in config_in.as_dict().items():
            self.config[section].update(vals)



    @staticmethod
    def print_version():
        """print doit version (includes path location)"""
        print(".".join([str(i) for i in VERSION]))
        print("lib @", os.path.dirname(os.path.abspath(__file__)))


    def get_cmds(self):
        """get all sub-commands
        :return dict: name - Command class
        """
        sub_cmds = PluginDict()
        # core doit commands
        for cmd_cls in self.DOIT_CMDS:
            sub_cmds[cmd_cls.get_name()] = cmd_cls
        # plugin commands
        sub_cmds.add_plugins(self.config, 'COMMAND')
        return sub_cmds


    def process_args(self, cmd_args):
        """process cmd line set "global" variables/parameters
        return list of args without processed variables
        """
        # get cmdline variables from args
        reset_vars()
        args_no_vars = []
        for arg in cmd_args:
            if (arg[0] != '-') and ('=' in arg):
                name, value = arg.split('=', 1)
                set_var(name, value)
            else:
                args_no_vars.append(arg)
        return args_no_vars


    def get_commands(self):  # pragma: no cover
        '''Notice for application subclassing DoitMain with old API'''
        msg = ('ERROR: You are using doit version {}, it is too new! '
               'This application requires version <= 0.27.\n')
        sys.stderr.write(msg.format('.'.join(str(v) for v in VERSION)))
        sys.exit(3)


    def run(self, all_args):
        """entry point for all commands

        :param all_args: list of string arguments from command line

        return codes:
          0: tasks executed successfully
          1: one or more tasks failed
          2: error while executing a task
          3: error before task execution starts,
             in this case the Reporter is not used.
             So be aware if you expect a different formatting (like JSON)
             from the Reporter.
        """
        # get list of available commands
        sub_cmds = self.get_cmds()
        task_loader = get_loader(self.config, self.task_loader, sub_cmds)

        # special parameters that dont run anything
        if all_args:
            if all_args[0] == "--version":
                self.print_version()
                return 0
            if all_args[0] == "--help":
                help = Help(task_loader=task_loader,
                            config=self.config,
                            bin_name=self.BIN_NAME,
                            cmds=sub_cmds)
                help.print_usage(sub_cmds.to_dict())
                return 0

        # loader command options might appear before command name
        try:
            loader_opt_parser = CmdParse(
                [CmdOption(opt) for opt in task_loader.cmd_options])
            loader_params, cmd_args = loader_opt_parser.parse_only(all_args)
        except CmdParseError:
            # normal to fail parsing if RUN command is not explicit
            loader_params = {}
            cmd_args = all_args

        # get "global vars" from cmd-line
        args = self.process_args(cmd_args)

        # get specified sub-command or use default='run'
        if len(args) == 0 or args[0] not in sub_cmds:
            specified_run = False
            cmd_name = 'run'
        else:
            specified_run = True
            cmd_name = args.pop(0)

        # execute command
        command = sub_cmds.get_plugin(cmd_name)(
            task_loader=task_loader,
            config=self.config,
            bin_name=self.BIN_NAME,
            cmds=sub_cmds,
            opt_vals=loader_params,
        )

        try:
            return command.parse_execute(args)

        # dont show traceback for user errors.
        except (CmdParseError, InvalidDodoFile,
                InvalidCommand, InvalidTask) as err:
            if isinstance(err, InvalidCommand):
                err.cmd_used = cmd_name if specified_run else None
                err.bin_name = self.BIN_NAME
            sys.stderr.write("ERROR: %s\n" % str(err))
            return 3

        except Exception:
            if command.pdb:  # pragma: no cover
                import pdb
                pdb.post_mortem(sys.exc_info()[2])
            sys.stderr.write(traceback.format_exc())
            return 3
