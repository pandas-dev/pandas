"""APIs to execute doit in non standard way.

- run(): shortcut to get tasks from current module/dict (instead of dodo.py).
- run_tasks(): to be used by custom CLIs, run tasks without CLI parsing.

"""

import sys

from doit.cmdparse import CmdParseError
from doit.exceptions import InvalidDodoFile, InvalidCommand, InvalidTask
from doit.cmd_base import ModuleTaskLoader, get_loader
from doit.doit_cmd import DoitMain


def run(task_creators):
    """run doit using task_creators

    @param task_creators: module or dict containing task creators
    """
    sys.exit(DoitMain(ModuleTaskLoader(task_creators)).run(sys.argv[1:]))




def run_tasks(loader, tasks, extra_config=None):
    """run DoitMain instance with speficied tasks and parameters

    :params tasks: list of task names (str)
    """
    loader.task_opts = tasks  # task_opts will be used as @task_param
    main = DoitMain(loader, extra_config=extra_config)
    task_names = list(tasks.keys())

    # get list of available commands
    sub_cmds = main.get_cmds()
    task_loader = get_loader(main.config, main.task_loader, sub_cmds)

    # execute command
    cmd_name = 'run'
    command = sub_cmds.get_plugin(cmd_name)(
        task_loader=task_loader,
        config=main.config,
        bin_name=main.BIN_NAME,
        cmds=sub_cmds,
        opt_vals={},
    )

    try:
        return command.parse_execute(task_names)
    except (CmdParseError, InvalidDodoFile,
            InvalidCommand, InvalidTask) as err:
        if isinstance(err, InvalidCommand):
            err.cmd_used = cmd_name
            err.bin_name = main.BIN_NAME
        raise err
