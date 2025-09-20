"""command doit info - display info on task metadata"""

import pprint

from .cmd_base import DoitCmdBase
from .exceptions import InvalidCommand


opt_hide_status = {
    'name': 'hide_status',
    'long': 'no-status',
    'type': bool,
    'default': False,
    'help': """Hides reasons why this task would be executed.
 [default: %(default)s]"""
}


class Info(DoitCmdBase):
    """command doit info"""

    doc_purpose = "show info about a task"
    doc_usage = "TASK"
    doc_description = None

    cmd_options = (opt_hide_status, )

    def _execute(self, pos_args, hide_status=False):
        if len(pos_args) != 1:
            msg = ('`info` failed, must select *one* task.'
                   '\nCheck `{} help info`.'.format(self.bin_name))
            raise InvalidCommand(msg)

        task_name = pos_args[0]
        # dict of all tasks
        tasks = dict([(t.name, t) for t in self.task_list])

        printer = pprint.PrettyPrinter(indent=4, stream=self.outstream)

        task = tasks[task_name]
        task_attrs = (
            ('file_dep', 'list'),
            ('task_dep', 'list'),
            ('setup_tasks', 'list'),
            ('calc_dep', 'list'),
            ('targets', 'list'),
            # these fields usually contains reference to python functions
            # 'actions', 'clean', 'uptodate', 'teardown', 'title'
            ('getargs', 'dict'),
            ('params', 'list'),
            ('verbosity', 'scalar'),
            ('watch', 'list'),
            ('meta', 'dict')
        )

        self.outstream.write('\n{}\n'.format(task.name))
        if task.doc:
            self.outstream.write('\n{}\n'.format(task.doc))

        # print reason task is not up-to-date
        retcode = 0
        if not hide_status:
            status = self.dep_manager.get_status(task, tasks, get_log=True)
            self.outstream.write('\n{:11s}: {}\n'
                                 .format('status', status.status))
            if status.status != 'up-to-date':
                # status.status == 'run' or status.status == 'error'
                self.outstream.write(self.get_reasons(status.reasons))
                self.outstream.write('\n')
                retcode = 1

        for (attr, attr_type) in task_attrs:
            value = getattr(task, attr)
            # only print fields that have non-empty value
            if value:
                self.outstream.write('\n{:11s}: '.format(attr))
                if attr_type == 'list':
                    self.outstream.write('\n')
                    for val in value:
                        self.outstream.write(' - {}\n'.format(val))
                else:
                    printer.pprint(getattr(task, attr))

        return retcode

    @staticmethod
    def get_reasons(reasons):
        '''return string with description of reason task is not up-to-date'''
        lines = []
        if reasons['has_no_dependencies']:
            lines.append(' * The task has no dependencies.')

        if reasons['uptodate_false']:
            lines.append(' * The following uptodate objects evaluate to false:')
            for utd, utd_args, utd_kwargs in reasons['uptodate_false']:
                msg = '    - {} (args={}, kwargs={})'
                lines.append(msg.format(utd, utd_args, utd_kwargs))

        if reasons['checker_changed']:
            msg = ' * The file_dep checker changed from {0} to {1}.'
            lines.append(msg.format(*reasons['checker_changed']))

        sentences = {
            'missing_target': 'The following targets do not exist:',
            'changed_file_dep': 'The following file dependencies have changed:',
            'missing_file_dep': 'The following file dependencies are missing:',
            'removed_file_dep': 'The following file dependencies were removed:',
            'added_file_dep': 'The following file dependencies were added:',
        }
        for reason, sentence in sentences.items():
            entries = reasons.get(reason)
            if entries:
                lines.append(' * {}'.format(sentence))
                for item in entries:
                    lines.append('    - {}'.format(item))
        return '\n'.join(lines)
