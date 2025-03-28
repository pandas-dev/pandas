import fnmatch
from collections import OrderedDict

from .control import TaskControl
from .cmd_base import DoitCmdBase
from .cmd_base import check_tasks_exist


opt_clean_dryrun = {
    'name': 'dryrun',
    'short': 'n',  # like make dry-run
    'long': 'dry-run',
    'type': bool,
    'default': False,
    'help': 'print actions without really executing them',
}

opt_clean_cleandep = {
    'name': 'cleandep',
    'short': 'c',
    'long': 'clean-dep',
    'type': bool,
    'default': False,
    'help': 'clean task dependencies too',
}

opt_clean_cleanall = {
    'name': 'cleanall',
    'short': 'a',  # all
    'long': 'clean-all',
    'type': bool,
    'default': False,
    'help': 'clean all task',
}

opt_clean_forget = {
    'name': 'cleanforget',
    'long': 'forget',
    'type': bool,
    'default': False,
    'help': 'also forget tasks after cleaning',
}


class Clean(DoitCmdBase):
    doc_purpose = "clean action / remove targets"
    doc_usage = "[TASK ...]"
    doc_description = ("If no task is specified clean default tasks and "
                       "set --clean-dep automatically.")

    cmd_options = (opt_clean_cleandep, opt_clean_cleanall,
                   opt_clean_dryrun, opt_clean_forget)


    def clean_tasks(self, tasks, dryrun, cleanforget):
        """ensure task clean-action is executed only once"""
        cleaned = set()
        forget_tasks = cleanforget and not dryrun
        for task in tasks:
            if task.name not in cleaned:
                cleaned.add(task.name)
                task.clean(self.outstream, dryrun)
                if forget_tasks:
                    self.dep_manager.remove(task.name)
        self.dep_manager.close()

    def _expand(self, clean_list):
        result = []
        for name in clean_list:
            if '*' in name:
                result.extend(t.name for t in self.task_list
                              if fnmatch.fnmatch(t.name, name))
            else:
                result.append(name)
        return result

    def _execute(self, dryrun, cleandep, cleanall, cleanforget,
                 pos_args=None):
        """Clean tasks
        @param task_list (list - L{Task}): list of all tasks from dodo file
        @ivar dryrun (bool): if True clean tasks are not executed
                            (just print out what would be executed)
        @param cleandep (bool): execute clean from task_dep
        @param cleanall (bool): clean all tasks
        @param cleanforget (bool): forget cleaned tasks
        @var default_tasks (list - string): list of default tasks
        @var selected_tasks (list - string): list of tasks selected
                                             from cmd-line
        """
        tasks = TaskControl(self.task_list).tasks
        # behavior of cleandep is different if selected_tasks comes from
        # command line or DOIT_CONFIG.default_tasks
        selected_tasks = pos_args
        check_tasks_exist(tasks, selected_tasks, skip_wildcard=True)

        # get base list of tasks to be cleaned
        if selected_tasks and not cleanall:  # from command line
            clean_list = self._expand(selected_tasks)
        else:
            # if not cleaning specific task enable clean_dep automatically
            cleandep = True
            if self.sel_tasks is not None:
                clean_list = self._expand(self.sel_tasks)  # default tasks from config
            else:
                clean_list = [t.name for t in self.task_list]
            # note: reversing is not required, but helps reversing
            # execution order even if there are no restrictions about order.
            clean_list.reverse()

        tree = CleanDepTree()
        # include dependencies in list
        if cleandep:
            for name in clean_list:
                tree.build_nodes_with_deps(tasks, name)
        # include only subtasks in list
        else:
            tree.build_nodes(tasks, clean_list)

        to_clean = [tasks[x] for x in tree.flat()]
        self.clean_tasks(to_clean, dryrun, cleanforget)


class CleanDepTree:
    """Create node structure where each node is a task and its children
    are tasks that has the node as a task_dep/setup_task.
    This creates an upside-down tree where leaf nodes should be
    the first ones to be "cleaned".
    """
    def __init__(self):
        self.nodes = OrderedDict()
        self._processed = set()  # task names that were already built

    def build_nodes_with_deps(self, tasks, task_name):
        """build node including task_dep's"""
        if task_name in self._processed:
            return
        else:
            self._processed.add(task_name)

        # add node itself if not in list of nodes
        self.nodes.setdefault(task_name, [])
        task = tasks[task_name]
        # reversing not required
        for dep_name in reversed(task.setup_tasks + task.task_dep):
            rev_dep = self.nodes.setdefault(dep_name, [])
            rev_dep.append(task_name)
            self.build_nodes_with_deps(tasks, dep_name)

    def build_nodes(self, tasks, clean_list):
        """build nodes with sub-tasks but no other task_dep"""
        for name in clean_list:
            # add node itself if not in list of nodes
            self.nodes.setdefault(name, [])
            task = tasks[name]
            # reversing not required
            for dep_name in reversed(task.task_dep):
                if tasks[dep_name].subtask_of == name:
                    rev_dep = self.nodes.setdefault(dep_name, [])
                    rev_dep.append(name)

    def flat(self):
        """return list of tasks in the order they should be `clean` """
        to_clean = []
        while self.nodes:
            head, children = self.nodes.popitem(0)
            to_clean.extend([x for x in self._get_leafs(head, children)])
        return to_clean

    def _get_leafs(self, name, children):
        for child_name in children:
            if child_name in self.nodes:
                grand = self.nodes.pop(child_name)
                yield from self._get_leafs(child_name, grand)
        yield name
