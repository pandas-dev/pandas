from .cmd_base import DoitCmdBase, check_tasks_exist, subtasks_iter


class Ignore(DoitCmdBase):
    doc_purpose = "ignore task (skip) on subsequent runs"
    doc_usage = "TASK [TASK ...]"
    doc_description = None

    cmd_options = ()

    def _execute(self, pos_args):
        """mark tasks to be ignored
        @param ignore_tasks: (list - str) tasks to be ignored.
        """
        ignore_tasks = pos_args
        # no task specified.
        if not ignore_tasks:
            msg = "You cant ignore all tasks! Please select a task.\n"
            self.outstream.write(msg)
            return

        tasks = dict([(t.name, t) for t in self.task_list])
        check_tasks_exist(tasks, ignore_tasks)

        for task_name in ignore_tasks:
            # for group tasks also remove all tasks from group
            sub_list = [t.name for t in subtasks_iter(tasks, tasks[task_name])]
            for to_ignore in [task_name] + sub_list:
                # ignore it - remove from dependency file
                self.dep_manager.ignore(tasks[to_ignore])
                self.outstream.write("ignoring %s\n" % to_ignore)

        self.dep_manager.close()
