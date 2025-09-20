from .cmd_base import DoitCmdBase, check_tasks_exist
from .cmd_base import subtasks_iter
import os


class ResetDep(DoitCmdBase):
    name = "reset-dep"
    doc_purpose = ("recompute and save the state of file dependencies without "
                   "executing actions")
    doc_usage = "[TASK ...]"
    cmd_options = ()
    doc_description = """
This command allows to recompute the information on file dependencies
(timestamp, md5sum, ... depending on the ``check_file_uptodate`` setting), and
save this in the database, without executing the actions.

The command run on all tasks by default, but it is possible to specify a list
of tasks to work on.

This is useful when the targets of your tasks already exist, and you want doit
to consider your tasks as up-to-date. One use-case for this command is when you
change the ``check_file_uptodate`` setting, which cause doit to consider all
your tasks as not up-to-date. It is also useful if you start using doit while
some of your data as already been computed, or when you add a file dependency
to a task that has already run.
"""

    def _execute(self, pos_args=None):
        filter_tasks = pos_args

        # dict of all tasks
        tasks = dict([(t.name, t) for t in self.task_list])

        # select tasks that command will be applied to
        if filter_tasks:
            # list only tasks passed on command line
            check_tasks_exist(tasks, filter_tasks)
            # get task by name
            task_list = []
            for name in filter_tasks:
                task = tasks[name]
                task_list.append(task)
                task_list.extend(subtasks_iter(tasks, task))
        else:
            task_list = self.task_list

        write = self.outstream.write
        for task in task_list:
            # Get these now because dep_manager.get_status will remove the task
            # from the db if the checker changed.
            values = self.dep_manager.get_values(task.name)
            result = self.dep_manager.get_result(task.name)

            missing_deps = [dep for dep in task.file_dep
                            if not os.path.exists(dep)]

            if len(missing_deps) > 0:
                deps = "', '".join(missing_deps)
                write(f"failed {task.name} (Dependent file '{deps}' does not exist.)\n")
                continue

            res = self.dep_manager.get_status(task, tasks)

            # An 'up-to-date' status means that it is useless to recompute the
            # state: file deps and targets exists, the state has not changed,
            # there is nothing more to do.
            if res.status == 'up-to-date':
                write("skip {}\n".format(task.name))
                continue

            task.values = values
            self.dep_manager.save_success(task, result_hash=result)
            write("processed {}\n".format(task.name))

        self.dep_manager.close()
