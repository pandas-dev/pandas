"""generate shell script with tab completion code for doit commands/tasks"""

import sys
from string import Template

from .exceptions import InvalidCommand
from .cmd_base import DoitCmdBase

opt_shell = {
    'name': 'shell',
    'short': 's',
    'long': 'shell',
    'type': str,
    'choices': (('bash', ''), ('zsh', '')),
    'default': 'bash',
    'help': 'Completion code for SHELL. [default: %(default)s]',
}

opt_hardcode_tasks = {
    'name': 'hardcode_tasks',
    'short': '',
    'long': 'hardcode-tasks',
    'type': bool,
    'default': False,
    'help': 'Hardcode tasks from current task list.',
}



class TabCompletion(DoitCmdBase):
    """generate scripts for tab-completion

    If hardcode-tasks options is chosen it will get the task
    list from the current dodo file and include in the completion script.
    Otherwise the script will dynamically call `doit list` to get the list
    of tasks.

    If it is completing a sub-task (contains ':' in the name),
    it will always call doit while evaluating the options.

    """
    doc_purpose = "generate script for tab-completion"
    doc_usage = ""
    doc_description = None

    cmd_options = (opt_shell, opt_hardcode_tasks, )

    def __init__(self, cmds=None, **kwargs):
        super(TabCompletion, self).__init__(cmds=cmds, **kwargs)
        self.init_kwargs = kwargs
        self.init_kwargs['cmds'] = cmds
        if cmds:
            self.cmds = cmds.to_dict()  # dict name - Command class

    def execute(self, opt_values, pos_args):
        if opt_values['shell'] == 'bash':
            self._generate_bash(opt_values, pos_args)
        elif opt_values['shell'] == 'zsh':
            self._generate_zsh(opt_values, pos_args)
        else:
            msg = 'Invalid option for --shell "{0}"'
            raise InvalidCommand(msg.format(opt_values['shell']))

    @classmethod
    def _bash_cmd_args(cls, cmd):
        """return case item for completion of specific sub-command"""
        comp = []
        if 'TASK' in cmd.doc_usage:
            comp.append('${tasks}')
        if 'COMMAND' in cmd.doc_usage:
            comp.append('${sub_cmds}')
        if comp:
            completion = '-W "{0}"'.format(' '.join(comp))
        else:
            completion = '-f'  # complete file
        return bash_subcmd_arg.format(cmd_name=cmd.name, completion=completion)


    def _generate_bash(self, opt_values, pos_args):
        # some applications built with doit do not use dodo.py files
        for opt in self.get_options():
            if opt.name == 'dodoFile':
                get_dodo_part = bash_get_dodo
                pt_list_param = '--file="$dodof"'
                break
        else:
            get_dodo_part = ''
            pt_list_param = ''

        # dict with template values
        pt_bin_name = sys.argv[0].split('/')[-1]
        tmpl_vars = {
            'pt_bin_name': pt_bin_name,
            'pt_cmds': ' '.join(sorted(self.cmds)),
            'pt_list_param': pt_list_param,
        }

        # if hardcode tasks
        if opt_values['hardcode_tasks']:
            if getattr(self.loader, 'API', 1) == 2:
                self.loader.setup(opt_values)
                self.loader.load_doit_config()
                self.task_list = self.loader.load_tasks(cmd=self, pos_args=pos_args)
            else:
                self.task_list, _ = self.loader.load_tasks(
                    self, opt_values, pos_args)
            task_names = (t.name for t in self.task_list if not t.subtask_of)
            tmpl_vars['pt_tasks'] = '"{0}"'.format(' '.join(sorted(task_names)))
        else:
            tmpl_list_cmd = "$({0} list {1} --quiet 2>/dev/null)"
            tmpl_vars['pt_tasks'] = tmpl_list_cmd.format(pt_bin_name,
                                                         pt_list_param)

        # case statement to complete sub-commands
        cmds_args = []
        for name in sorted(self.cmds):
            cmd_class = self.cmds[name]
            cmd = cmd_class(**self.init_kwargs)
            cmds_args.append(self._bash_cmd_args(cmd))
        comp_subcmds = ("\n    case ${words[1]} in\n"
                        + "".join(cmds_args)
                        + "\n    esac\n")

        template = Template(
            bash_start + bash_opt_file + get_dodo_part
            + bash_task_list + bash_first_arg
            + comp_subcmds + bash_end)
        self.outstream.write(template.safe_substitute(tmpl_vars))


    @staticmethod
    def _zsh_arg_line(opt):
        """create a text line for completion of a command arg"""
        # '(-c|--continue)'{-c,--continue}'[continue executing tasks...]' \
        # '--db-file[file used to save successful runs]' \
        if opt.short and opt.long:
            tmpl = ('"(-{0.short}|--{0.long})"{{-{0.short},--{0.long}}}"'
                    '[{help}]" \\')
        elif not opt.short and opt.long:
            tmpl = '"--{0.long}[{help}]" \\'
        elif opt.short and not opt.long:
            tmpl = '"-{0.short}[{help}]" \\'
        else:  # without short or long options cant be really used
            return ''
        ohelp = opt.help.replace(']', r'\]').replace('"', r'\"')
        return tmpl.format(opt, help=ohelp).replace('\n', ' ')


    @classmethod
    def _zsh_arg_list(cls, cmd):
        """return list of arguments lines for zsh completion"""
        args = []
        for opt in cmd.get_options():
            args.append(cls._zsh_arg_line(opt))
        if 'TASK' in cmd.doc_usage:
            args.append("'*::task:(($tasks))'")
        if 'COMMAND' in cmd.doc_usage:
            args.append("'::cmd:(($commands))'")
        return args

    @classmethod
    def _zsh_cmd_args(cls, cmd):
        """create the content for "case" statement with all command options """
        arg_lines = cls._zsh_arg_list(cmd)
        tmpl = """
      ({cmd_name})
          _command_args=(
            {args_body}
            ''
        )
      ;;
"""
        args_body = '\n            '.join(arg_lines)
        return tmpl.format(cmd_name=cmd.name, args_body=args_body)


    # TODO:
    # detect correct dodo-file location
    # complete sub-tasks
    # task options
    def _generate_zsh(self, opt_values, pos_args):
        # deal with doit commands
        cmds_desc = []
        cmds_args = []
        for name in sorted(self.cmds):
            cmd_class = self.cmds[name]
            cmd = cmd_class(**self.init_kwargs)
            cmds_desc.append("    '{0}: {1}'".format(cmd.name, cmd.doc_purpose))
            cmds_args.append(self._zsh_cmd_args(cmd))

        template_vars = {
            'pt_bin_name': sys.argv[0].split('/')[-1],
            'pt_cmds': '\n    '.join(cmds_desc),
            'pt_cmds_args': '\n'.join(cmds_args),
        }

        if opt_values['hardcode_tasks']:
            if getattr(self.loader, 'API', 1) == 2:
                self.loader.setup(opt_values)
                self.loader.load_doit_config()
                self.task_list = self.loader.load_tasks(cmd=self, pos_args=pos_args)
            else:
                self.task_list, _ = self.loader.load_tasks(
                    self, opt_values, pos_args)
            lines = []
            for task in self.task_list:
                if not task.subtask_of:
                    lines.append("'{0}: {1}'".format(task.name, task.doc))
            template_vars['pt_tasks'] = '(\n{0}\n)'.format(
                '\n'.join(sorted(lines)))
        else:
            tmp_tasks = Template(
                '''("${(f)$($pt_bin_name list --template '{name}: {doc}')}")''')
            template_vars['pt_tasks'] = tmp_tasks.safe_substitute(template_vars)


        template = Template(zsh_start)
        self.outstream.write(template.safe_substitute(template_vars))




############## templates
# Variables starting with 'pt_' belongs to the Python Template
# to generate the script.
# Remaining are shell variables used in the script.


################################################################
############### bash template


bash_start = """# bash completion for $pt_bin_name
# auto-generate by `$pt_bin_name tabcompletion`

# to activate it you need to 'source' the generate script
# $ source <generated-script>

# reference => http://www.debian-administration.org/articles/317
# patch => http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=711879

_$pt_bin_name()
{
    local cur prev words cword basetask sub_cmds tasks i dodof
    COMPREPLY=() # contains list of words with suitable completion
    # remove colon from word separator list because doit uses colon on task names
    _get_comp_words_by_ref -n : cur prev words cword
    # list of sub-commands
    sub_cmds="$pt_cmds"

"""

# FIXME - wont be necessary after adding support for options with type
bash_opt_file = """
    # options that take file/dir as values should complete file-system
    if [[ "$prev" == "-f" || "$prev" == "-d" || "$prev" == "-o" ]]; then
        _filedir
        return 0
    fi
    if [[ "$cur" == *=* ]]; then
        prev=${cur/=*/}
        cur=${cur/*=/}
        if [[ "$prev" == "--file=" || "$prev" == "--dir=" || "$prev" == "--output-file=" ]]; then
            _filedir -o nospace
            return 0
        fi
    fi

"""


bash_get_dodo = """
    # get name of the dodo file
    for (( i=0; i < ${#words[@]}; i++)); do
        case "${words[i]}" in
        -f)
            dodof=${words[i+1]}
            break
            ;;
        --file=*)
            dodof=${words[i]/*=/}
            break
            ;;
        esac
    done
    # dodo file not specified, use default
    if [ ! $dodof ]
      then
         dodof="dodo.py"
    fi

"""

bash_task_list = """
    # get task list
    # if it there is colon it is getting a subtask, complete only subtask names
    if [[ "$cur" == *:* ]]; then
        # extract base task name (remove everything after colon)
        basetask=${cur%:*}
        # sub-tasks
        tasks=$($pt_bin_name list $pt_list_param --quiet --all ${basetask} 2>/dev/null)
        COMPREPLY=( $(compgen -W "${tasks}" -- ${cur}) )
        __ltrim_colon_completions "$cur"
        return 0
    # without colons get only top tasks
    else
        tasks=$pt_tasks
    fi

"""

bash_first_arg = """
    # match for first parameter must be sub-command or task
    # FIXME doit accepts options "-" in the first parameter but we ignore this case
    if [[ ${cword} == 1 ]] ; then
        COMPREPLY=( $(compgen -W "${sub_cmds} ${tasks}" -- ${cur}) )
        return 0
    fi
"""

bash_subcmd_arg = """
        {cmd_name})
            COMPREPLY=( $(compgen {completion} -- $cur) )
            return 0
            ;;"""

bash_end = """
    # if there is already one parameter match only tasks (no commands)
    COMPREPLY=( $(compgen -W "${tasks}" -- ${cur}) )

}
complete -o filenames -F _$pt_bin_name $pt_bin_name
"""



################################################################
############### zsh template


zsh_start = """#compdef $pt_bin_name

_$pt_bin_name() {
    local -a commands tasks
    # format is 'completion:description'
    commands=(
    $pt_cmds
    )

    # split output by lines to create an array
    tasks=$pt_tasks

    # complete command or task name
    if (( CURRENT == 2 )); then
        _arguments -A : '::cmd:(($commands))' '::task:(($tasks))'
        return
    fi

    # revome program name from $words and decrement CURRENT
    local curcontext context state state_desc line
    _arguments -C '*:: :->'

    # complete sub-command or task options
    local -a _command_args
    case "$words[1]" in
        $pt_cmds_args

        # default completes task names
        (*)
           _command_args='*::task:(($tasks))'
        ;;
    esac

    # -A no options will be completed after the first non-option argument
    _arguments -A : $_command_args
    return 0
}

_$pt_bin_name
"""
