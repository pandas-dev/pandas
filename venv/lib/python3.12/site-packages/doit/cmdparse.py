"""Parse command line options and execute it.

Built on top of getopt. optparse can't handle sub-commands.
"""
import os
import getopt
import copy
from collections import OrderedDict



class DefaultUpdate(dict):
    """A dictionary that has an "update_defaults" method where
    only items with default values are updated.

    This is used when you have a dict that has multiple source of values
    (i.e. hardcoded, config file, command line). And values are updated
    beginning from the source with higher priority.

    A default value is added with the method set_default or add_defaults.
    """
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        # set of keys that have a non-default value
        self._non_default_keys = set()

    def set_default(self, key, value):
        """set default value for given key"""
        dict.__setitem__(self, key, value)

    def add_defaults(self, source):
        """add default values from another dict
        @param source: (dict)"""
        for key, value in source.items():
            if key not in self:
                self.set_default(key, value)

    def update_defaults(self, update_dict):
        """like dict.update but do not update items that have
        a non-default value"""
        for key, value in update_dict.items():
            if key in self._non_default_keys:
                continue
            self.set_default(key, value)

    def __setitem__(self, key, value):
        """overwrite to keep track of _non_default_keys"""
        try:
            self._non_default_keys.add(key)
        # http://bugs.python.org/issue826897
        except AttributeError:
            self._non_default_keys = set()
            self._non_default_keys.add(key)
        dict.__setitem__(self, key, value)


class CmdParseError(Exception):
    """Error parsing options """


class CmdOption(object):
    """a command line option

       - name (string) : variable name
       - section (string): meta info used to group entries when generating help
       - default (value from its type): default value
       - type (type): type of the variable. must be able to be initialized
                      taking a single string parameter.
                      if type is bool. option is just a flag. and if present
                      its value is set to True.
       - short (string): argument short name
       - long (string): argument long name
       - inverse (string): argument long name to be the inverse of the default
                           value (only used by boolean options)
       - choices(list - 2-tuple str): sequence of 2-tuple of choice name,
                                      choice description.
       - help (string): option description
    """

    def __init__(self, opt_dict):
        # options must contain 'name' and 'default' value
        opt_dict = opt_dict.copy()
        for field in ('name', 'default',):
            if field not in opt_dict:
                msg = "CmdOption dict %r missing required property '%s'"
                raise CmdParseError(msg % (opt_dict, field))

        self.name = opt_dict.pop('name')
        self.section = opt_dict.pop('section', '')
        self.type = opt_dict.pop('type', str)
        self.set_default(opt_dict.pop('default'))
        self.short = opt_dict.pop('short', '')
        self.long = opt_dict.pop('long', '')
        self.inverse = opt_dict.pop('inverse', '')
        self.choices = dict(opt_dict.pop('choices', []))
        self.help = opt_dict.pop('help', '')
        self.metavar = opt_dict.pop('metavar', 'ARG')
        self.env_var = opt_dict.pop('env_var', None)

        # TODO add some hint for tab-completion scripts

        # options can not contain any unrecognized field
        if opt_dict:
            msg = "CmdOption dict contains invalid property '%s'"
            raise CmdParseError(msg % list(opt_dict.keys()))

    def __repr__(self):
        tmpl = ("{0}({{'name':{1.name!r}, "
                "'short':{1.short!r},"
                "'long':{1.long!r} }})")
        return tmpl.format(self.__class__.__name__, self)

    def set_default(self, val):
        """set default value, value is already the expected type"""
        if self.type is list:
            self.default = copy.copy(val)
        else:
            self.default = val

    def validate_choice(self, given_value):
        """raise error is value is not a valid choice"""
        if given_value not in self.choices:
            msg = ("Error parsing parameter '{}'. "
                   "Provided '{}' but available choices are: {}.")
            choices = ", ".join(f"'{k}'" for k in self.choices.keys())
            raise CmdParseError(msg.format(self.name, given_value, choices))


    _boolean_states = {
        '1': True, 'yes': True, 'true': True, 'on': True,
        '0': False, 'no': False, 'false': False, 'off': False,
    }
    def str2boolean(self, str_val):
        """convert string to boolean"""
        try:
            return self._boolean_states[str_val.lower()]
        except Exception:
            raise ValueError('Not a boolean: {}'.format(str_val))

    def str2type(self, str_val):
        """convert string value to option type value"""
        try:
            # no conversion if value is not a string
            if not isinstance(str_val, str):
                val = str_val
            elif self.type is bool:
                val = self.str2boolean(str_val)
            elif self.type is list:
                parts = [p.strip() for p in str_val.split(',')]
                val = [p for p in parts if p]  # remove empty strings
            else:
                val = self.type(str_val)
        except ValueError as exception:
            msg = (f"Error parsing parameter '{self.name}' {self.type}.\n"
                   f"{exception}\n")
            raise CmdParseError(msg)

        if self.choices:
            self.validate_choice(val)
        return val


    @staticmethod
    def _print_2_columns(col1, col2):
        """print using a 2-columns format """
        column1_len = 24
        column2_start = 28
        left = (col1).ljust(column1_len)
        right = col2.replace('\n', '\n' + column2_start * ' ')
        return "  %s  %s" % (left, right)

    def help_param(self):
        """return string of option's short and long name
        i.e.:   -f ARG, --file=ARG
        """
        opts_str = []
        if self.short:
            if self.type is bool:
                opts_str.append('-%s' % self.short)
            else:
                opts_str.append('-%s %s' % (self.short, self.metavar))
        if self.long:
            if self.type is bool:
                opts_str.append('--%s' % self.long)
            else:
                opts_str.append('--%s=%s' % (self.long, self.metavar))
        return ', '.join(opts_str)

    def help_choices(self):
        """return string with help for option choices"""
        if not self.choices:
            return ''

        # if choice has a description display one choice per line...
        if any(self.choices.values()):
            items = []
            for choice in sorted(self.choices):
                items.append("\n{}: {}".format(choice, self.choices[choice]))
            return "\nchoices:" + "".join(items)
        # ... otherwise display in a single line
        else:
            return "\nchoices: " + ", ".join(sorted(self.choices.keys()))


    def help_doc(self):
        """return list of string of option's help doc

        Note this is used only to display help on tasks.
        For commands a better and more complete version is used.
        see cmd_base:Command.help
        """
        # ignore option that cant be modified on cmd line
        if not (self.short or self.long):
            return []

        text = []
        opt_str = self.help_param()
        # TODO It should always display option's default value
        opt_help = self.help % {'default': self.default}
        opt_choices = self.help_choices()
        opt_config = 'config: {}'.format(self.name)
        opt_env = ', environ: {}'.format(self.env_var) if self.env_var else ''

        desc = f'{opt_help} {opt_choices} ({opt_config}{opt_env})'
        text.append(self._print_2_columns(opt_str, desc))
        # print bool inverse option
        if self.inverse:
            opt_str = '--%s' % self.inverse
            opt_help = 'opposite of --%s' % self.long
            text.append(self._print_2_columns(opt_str, opt_help))
        return text



class CmdParse(object):
    """Process string with command options

    @ivar options: (list - CmdOption)
    """
    _type = "Command"

    def __init__(self, options):
        self._options = OrderedDict((o.name, o) for o in options)

    def __contains__(self, key):
        return key in self._options

    def __getitem__(self, key):
        return self._options[key]

    @property
    def options(self):
        """return list of options for backward compatibility"""
        return list(self._options.values())

    def get_short(self):
        """return string with short options for getopt"""
        short_list = ""
        for opt in self._options.values():
            if not opt.short:
                continue
            short_list += opt.short
            # ':' means option takes a value
            if opt.type is not bool:
                short_list += ':'
        return short_list

    def get_long(self):
        """return list with long options for getopt"""
        long_list = []
        for opt in self._options.values():
            long_name = opt.long
            if not long_name:
                continue
            # '=' means option takes a value
            if opt.type is not bool:
                long_name += '='
            long_list.append(long_name)
            if opt.inverse:
                long_list.append(opt.inverse)
        return long_list

    def get_option(self, opt_str):
        """return tuple
            - CmdOption from matching opt_str. or None
            - (bool) matched inverse
        """
        for opt in self._options.values():
            if opt_str in ('-' + opt.short, '--' + opt.long):
                return opt, False
            if opt_str == '--' + opt.inverse:
                return opt, True
        return None, None

    def overwrite_defaults(self, new_defaults):
        """overwrite self.options default values

        This values typically come from an INI file
        """
        for key, val in new_defaults.items():
            if key in self._options:
                opt = self._options[key]
                opt.set_default(opt.str2type(val))


    def parse_only(self, in_args, params=None):
        """parse arguments into options(params) and positional arguments

        @param in_args (list - string): typically sys.argv[1:]
        @return params, args
             params(dict): params contain the actual values from the options.
                           where the key is the name of the option.
             pos_args (list - string): positional arguments
        """
        params = params if params else {}

        # parse cmdline options using getopt
        try:
            opts, args = getopt.getopt(in_args, self.get_short(),
                                       self.get_long())
        except Exception as error:
            msg = (f"Error parsing {self._type}: {error} "
                   f"(parsing options: {self.options}). Got: {in_args}")
            raise CmdParseError(msg)

        # update params with values from command line
        for opt, val in opts:
            this, inverse = self.get_option(opt)
            if this.type is bool:
                params[this.name] = not inverse
            elif this.type is list:
                params[this.name].append(val)
            else:
                params[this.name] = this.str2type(val)

        return params, args


    def parse(self, in_args):
        """parse arguments into options(params) and positional arguments

        Also get values from shell ENV.

        Returned params is a `DefaultUpdate` type and includes
        an item for every option.

        @param in_args (list - string): typically sys.argv[1:]
        @return params, args
             params(dict): params contain the actual values from the options.
                           where the key is the name of the option.
             pos_args (list - string): positional arguments
        """
        params = DefaultUpdate()
        # add default values
        for opt in self._options.values():
            params.set_default(opt.name, opt.default)

        # get values from shell ENV
        for opt in self._options.values():
            if opt.env_var:
                val = os.getenv(opt.env_var)
                if val is not None:
                    params[opt.name] = opt.str2type(val)

        return self.parse_only(in_args, params)



class TaskParse(CmdParse):
    """Process string with command options (for tasks)"""
    _type = "Task"
