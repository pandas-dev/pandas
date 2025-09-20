from __future__ import annotations

import functools
import inspect
import os
import sys
from collections import defaultdict
from textwrap import TextWrapper

from plumbum import colors, local
from plumbum.cli.i18n import get_translation_for
from plumbum.lib import getdoc

from .switches import (
    CountOf,
    Flag,
    MissingArgument,
    MissingMandatorySwitch,
    PositionalArgumentsError,
    SubcommandError,
    SwitchCombinationError,
    SwitchError,
    UnknownSwitch,
    WrongArgumentType,
    switch,
)
from .terminal import get_terminal_size

_translation = get_translation_for(__name__)
T_, ngettext = _translation.gettext, _translation.ngettext


class ShowHelp(SwitchError):
    pass


class ShowHelpAll(SwitchError):
    pass


class ShowVersion(SwitchError):
    pass


class SwitchParseInfo:
    __slots__ = ["swname", "val", "index", "__weakref__"]

    def __init__(self, swname, val, index):
        self.swname = swname
        self.val = val
        self.index = index


class Subcommand:
    def __init__(self, name, subapplication):
        self.name = name
        self.subapplication = subapplication

    def get(self):
        if isinstance(self.subapplication, str):
            modname, clsname = self.subapplication.rsplit(".", 1)
            mod = __import__(modname, None, None, "*")
            try:
                cls = getattr(mod, clsname)
            except AttributeError:
                raise ImportError(f"cannot import name {clsname}") from None
            self.subapplication = cls
        return self.subapplication

    def __repr__(self):
        return T_("Subcommand({self.name}, {self.subapplication})").format(self=self)


_switch_groups = ["Switches", "Meta-switches"]
_switch_groups_l10n = [T_("Switches"), T_("Meta-switches")]


# ===================================================================================================
# CLI Application base class
# ===================================================================================================


class Application:
    """The base class for CLI applications; your "entry point" class should derive from it,
    define the relevant switch functions and attributes, and the ``main()`` function.
    The class defines two overridable "meta switches" for version (``-v``, ``--version``)
    help (``-h``, ``--help``), and help-all (``--help-all``).

    The signature of the main function matters: any positional arguments (e.g., non-switch
    arguments) given on the command line are passed to the ``main()`` function; if you wish
    to allow unlimited number of positional arguments, use varargs (``*args``). The names
    of the arguments will be shown in the help message.

    The classmethod ``run`` serves as the entry point of the class. It parses the command-line
    arguments, invokes switch functions and enters ``main``. You should **not override** this
    method.

    Usage::

        class FileCopier(Application):
            stat = Flag("p", "copy stat info as well")

            def main(self, src, dst):
                if self.stat:
                    shutil.copy2(src, dst)
                else:
                    shutil.copy(src, dst)

        if __name__ == "__main__":
            FileCopier.run()

    There are several class-level attributes you may set:

    * ``PROGNAME`` - the name of the program; if ``None`` (the default), it is set to the
      name of the executable (``argv[0]``); can be in color. If only a color, will be applied to the name.

    * ``VERSION`` - the program's version (defaults to ``1.0``, can be in color)

    * ``DESCRIPTION`` - a short description of your program (shown in help). If not set,
      the class' ``__doc__`` will be used. Can be in color.

    * ``DESCRIPTION_MORE`` - a detailed description of your program (shown in help). The text will be printed
      by paragraphs (specified by empty lines between them). The indentation of each paragraph will be the
      indentation of its first line. List items are identified by their first non-whitespace character being
      one of '-', '*', and '/'; so that they are not combined with preceding paragraphs. Bullet '/' is
      "invisible", meaning that the bullet itself will not be printed to the output.

    * ``USAGE`` - the usage line (shown in help).

    * ``COLOR_USAGE_TITLE`` - The color of the usage line's header.

    * ``COLOR_USAGE`` - The color of the usage line.

    * ``COLOR_GROUPS`` - A dictionary that sets colors for the groups, like Meta-switches, Switches,
      and Subcommands.

    * ``COLOR_GROUP_TITLES`` - A dictionary that sets colors for the group titles. If the dictionary is empty,
      it defaults to ``COLOR_GROUPS``.

    * ``SUBCOMMAND_HELPMSG`` - Controls the printing of extra "see subcommand -h" help message.
      Default is a message, set to False to remove.

    * ``ALLOW_ABBREV`` - Controls whether partial switch names are supported, for example '--ver' will match
      '--verbose'. Default is False for backward consistency with previous plumbum releases. Note that ambiguous
      abbreviations will not match, for example if --foothis and --foothat are defined, then --foo will not match.

    A note on sub-commands: when an application is the root, its ``parent`` attribute is set to
    ``None``. When it is used as a nested-command, ``parent`` will point to its direct ancestor.
    Likewise, when an application is invoked with a sub-command, its ``nested_command`` attribute
    will hold the chosen sub-application and its command-line arguments (a tuple); otherwise, it
    will be set to ``None``

    """

    PROGNAME = None
    DESCRIPTION = None
    DESCRIPTION_MORE = None
    VERSION = None
    USAGE = None
    COLOR_USAGE = None
    COLOR_USAGE_TITLE = None
    COLOR_GROUPS = None
    COLOR_GROUP_TITLES = None
    CALL_MAIN_IF_NESTED_COMMAND = True
    SUBCOMMAND_HELPMSG = T_("see '{parent} {sub} --help' for more info")
    ALLOW_ABBREV = False

    parent = None
    nested_command = None
    _unbound_switches = ()

    def __new__(cls, executable=None):
        """Allows running the class directly as a shortcut for main.
        This is necessary for some setup scripts that want a single function,
        instead of an expression with a dot in it."""

        if executable is None:
            # This return value was not a class instance, so __init__ is never called
            return cls.run()

        return super().__new__(cls)

    def __init__(self, executable):
        # Filter colors

        if self.PROGNAME is None:
            self.PROGNAME = os.path.basename(executable)
        elif isinstance(self.PROGNAME, colors._style):
            self.PROGNAME = self.PROGNAME | os.path.basename(executable)
        elif colors.filter(self.PROGNAME) == "":
            self.PROGNAME = colors.extract(self.PROGNAME) | os.path.basename(executable)
        if self.DESCRIPTION is None:
            self.DESCRIPTION = getdoc(self)

        # Allow None for the colors
        self.COLOR_GROUPS = defaultdict(
            lambda: colors.do_nothing,
            {} if type(self).COLOR_GROUPS is None else type(self).COLOR_GROUPS,
        )

        self.COLOR_GROUP_TITLES = defaultdict(
            lambda: colors.do_nothing,
            self.COLOR_GROUPS
            if type(self).COLOR_GROUP_TITLES is None
            else type(self).COLOR_GROUP_TITLES,
        )
        if type(self).COLOR_USAGE is None:
            self.COLOR_USAGE = colors.do_nothing

        self.executable = executable
        self._switches_by_name = {}
        self._switches_by_func = {}
        self._switches_by_envar = {}
        self._subcommands = {}

        for cls in reversed(type(self).mro()):
            for obj in cls.__dict__.values():
                if isinstance(obj, Subcommand):
                    name = colors.filter(obj.name)
                    if name.startswith("-"):
                        raise SubcommandError(
                            T_("Sub-command names cannot start with '-'")
                        )
                    # it's okay for child classes to override sub-commands set by their parents
                    self._subcommands[name] = obj
                    continue

                swinfo = getattr(obj, "_switch_info", None)
                if not swinfo:
                    continue
                for name in swinfo.names:
                    if name in self._unbound_switches:
                        continue
                    if (
                        name in self._switches_by_name
                        and not self._switches_by_name[name].overridable
                    ):
                        raise SwitchError(
                            T_(
                                "Switch {name} already defined and is not overridable"
                            ).format(name=name)
                        )
                    self._switches_by_name[name] = swinfo
                    self._switches_by_func[swinfo.func] = swinfo
                    if swinfo.envname:
                        self._switches_by_envar[swinfo.envname] = swinfo

    @property
    def root_app(self):
        return self.parent.root_app if self.parent else self

    @classmethod
    def unbind_switches(cls, *switch_names):
        """Unbinds the given switch names from this application. For example

        ::

            class MyApp(cli.Application):
                pass
            MyApp.unbind_switches("--version")

        """
        cls._unbound_switches += tuple(
            name.lstrip("-") for name in switch_names if name
        )

    @classmethod
    def subcommand(cls, name, subapp=None):
        """Registers the given sub-application as a sub-command of this one. This method can be
        used both as a decorator and as a normal ``classmethod``::

            @MyApp.subcommand("foo")
            class FooApp(cli.Application):
                pass

        Or ::

            MyApp.subcommand("foo", FooApp)

        .. versionadded:: 1.1

        .. versionadded:: 1.3
            The sub-command can also be a string, in which case it is treated as a
            fully-qualified class name and is imported on demand. For example,

            MyApp.subcommand("foo", "fully.qualified.package.FooApp")

        """

        def wrapper(subapp):
            subname = subapp if isinstance(subapp, str) else subapp.__name__
            attrname = f"_subcommand_{subname}"
            setattr(cls, attrname, Subcommand(name, subapp))
            return subapp

        return wrapper(subapp) if subapp else wrapper

    def _get_partial_matches(self, partialname):
        matches = []
        for switch_ in self._switches_by_name:
            if switch_.startswith(partialname):
                matches += [
                    switch_,
                ]
        return matches

    def _parse_args(self, argv):
        tailargs = []
        swfuncs = {}
        index = 0

        while argv:
            index += 1
            a = argv.pop(0)
            val = None
            if a == "--":
                # end of options, treat the rest as tailargs
                tailargs.extend(argv)
                break

            if a in self._subcommands:
                subcmd = self._subcommands[a].get()
                self.nested_command = (
                    subcmd,
                    [self.PROGNAME + " " + self._subcommands[a].name, *argv],
                )
                break

            if a.startswith("--") and len(a) >= 3:
                # [--name], [--name=XXX], [--name, XXX], [--name, ==, XXX],
                # [--name=, XXX], [--name, =XXX]
                eqsign = a.find("=")
                if eqsign >= 0:
                    name = a[2:eqsign]
                    argv.insert(0, a[eqsign:])
                else:
                    name = a[2:]

                if self.ALLOW_ABBREV:
                    partials = self._get_partial_matches(name)
                    if len(partials) == 1:
                        name = partials[0]
                    elif len(partials) > 1:
                        raise UnknownSwitch(
                            T_("Ambiguous partial switch {0}").format("--" + name)
                        )

                swname = "--" + name
                if name not in self._switches_by_name:
                    raise UnknownSwitch(T_("Unknown switch {0}").format(swname))
                swinfo = self._switches_by_name[name]
                if swinfo.argtype:
                    if not argv:
                        raise MissingArgument(
                            T_("Switch {0} requires an argument").format(swname)
                        )
                    a = argv.pop(0)
                    if a and a[0] == "=":
                        if len(a) >= 2:
                            val = a[1:]
                        else:
                            if not argv:
                                raise MissingArgument(
                                    T_("Switch {0} requires an argument").format(swname)
                                )
                            val = argv.pop(0)
                    else:
                        val = a

            elif a.startswith("-") and len(a) >= 2:
                # [-a], [-a, XXX], [-aXXX], [-abc]
                name = a[1]
                swname = "-" + name
                if name not in self._switches_by_name:
                    raise UnknownSwitch(T_("Unknown switch {0}").format(swname))
                swinfo = self._switches_by_name[name]
                if swinfo.argtype:
                    if len(a) >= 3:
                        val = a[2:]
                    else:
                        if not argv:
                            raise MissingArgument(
                                T_("Switch {0} requires an argument").format(swname)
                            )
                        val = argv.pop(0)
                elif len(a) >= 3:
                    argv.insert(0, "-" + a[2:])

            else:
                if a.startswith("-"):
                    raise UnknownSwitch(T_("Unknown switch {0}").format(a))
                tailargs.append(a)
                continue

            # handle argument
            val = self._handle_argument(val, swinfo.argtype, name)

            if swinfo.func in swfuncs:
                if swinfo.list:
                    swfuncs[swinfo.func].val[0].append(val)
                else:
                    if swfuncs[swinfo.func].swname == swname:
                        raise SwitchError(T_("Switch {0} already given").format(swname))
                    raise SwitchError(
                        T_("Switch {0} already given ({1} is equivalent)").format(
                            swfuncs[swinfo.func].swname, swname
                        )
                    )
            else:
                if swinfo.list:
                    swfuncs[swinfo.func] = SwitchParseInfo(swname, ([val],), index)
                elif val is NotImplemented:
                    swfuncs[swinfo.func] = SwitchParseInfo(swname, (), index)
                else:
                    swfuncs[swinfo.func] = SwitchParseInfo(swname, (val,), index)

        # Extracting arguments from environment variables
        envindex = 0
        for env, swinfo in self._switches_by_envar.items():
            envindex -= 1
            envval = local.env.get(env)
            if envval is None:
                continue

            if swinfo.func in swfuncs:
                continue  # skip if overridden by command line arguments

            val = self._handle_argument(envval, swinfo.argtype, env)
            envname = f"${env}"
            if swinfo.list:
                # multiple values over environment variables are not supported,
                # this will require some sort of escaping and separator convention
                swfuncs[swinfo.func] = SwitchParseInfo(envname, ([val],), envindex)
            elif val is NotImplemented:
                swfuncs[swinfo.func] = SwitchParseInfo(envname, (), envindex)
            else:
                swfuncs[swinfo.func] = SwitchParseInfo(envname, (val,), envindex)

        return swfuncs, tailargs

    @classmethod
    def autocomplete(cls, argv):
        """This is supplied to make subclassing and testing argument completion methods easier"""

    @staticmethod
    def _handle_argument(val, argtype, name):
        if argtype:
            try:
                return argtype(val)
            except (TypeError, ValueError) as ex:
                raise WrongArgumentType(
                    T_(
                        "Argument of {name} expected to be {argtype}, not {val!r}:\n    {ex!r}"
                    ).format(name=name, argtype=argtype, val=val, ex=ex)
                ) from None
        else:
            return NotImplemented

    def _validate_args(self, swfuncs, tailargs):
        if self.help.__func__ in swfuncs:
            raise ShowHelp()
        if self.helpall.__func__ in swfuncs:
            raise ShowHelpAll()
        if self.version.__func__ in swfuncs:
            raise ShowVersion()

        requirements = {}
        exclusions = {}
        for swinfo in self._switches_by_func.values():
            if swinfo.mandatory and swinfo.func not in swfuncs:
                raise MissingMandatorySwitch(
                    T_("Switch {0} is mandatory").format(
                        "/".join(
                            ("-" if len(n) == 1 else "--") + n for n in swinfo.names
                        )
                    )
                )
            requirements[swinfo.func] = {
                self._switches_by_name[req] for req in swinfo.requires
            }
            exclusions[swinfo.func] = {
                self._switches_by_name[exc] for exc in swinfo.excludes
            }

        # TODO: compute topological order

        gotten = set(swfuncs.keys())
        for func in gotten:
            missing = {f.func for f in requirements[func]} - gotten
            if missing:
                raise SwitchCombinationError(
                    T_("Given {0}, the following are missing {1}").format(
                        swfuncs[func].swname,
                        [self._switches_by_func[f].names[0] for f in missing],
                    )
                )
            invalid = {f.func for f in exclusions[func]} & gotten
            if invalid:
                raise SwitchCombinationError(
                    T_("Given {0}, the following are invalid {1}").format(
                        swfuncs[func].swname, [swfuncs[f].swname for f in invalid]
                    )
                )

        m = inspect.getfullargspec(self.main)

        if sys.version_info < (3, 10):
            sig = inspect.signature(self.main)
        else:
            sig = inspect.signature(self.main, eval_str=True)

        max_args = sys.maxsize if m.varargs else len(m.args) - 1
        min_args = len(m.args) - 1 - (len(m.defaults) if m.defaults else 0)
        if len(tailargs) < min_args:
            raise PositionalArgumentsError(
                ngettext(
                    "Expected at least {0} positional argument, got {1}",
                    "Expected at least {0} positional arguments, got {1}",
                    min_args,
                ).format(min_args, tailargs)
            )
        if len(tailargs) > max_args:
            raise PositionalArgumentsError(
                ngettext(
                    "Expected at most {0} positional argument, got {1}",
                    "Expected at most {0} positional arguments, got {1}",
                    max_args,
                ).format(max_args, tailargs)
            )

        # Positional argument validation
        if hasattr(self.main, "positional"):
            tailargs = self._positional_validate(
                tailargs,
                self.main.positional,
                self.main.positional_varargs,
                m.args[1:],
                m.varargs,
            )

        elif hasattr(m, "annotations") and m.annotations:
            args_names = list(m.args[1:])
            positional = [None] * len(args_names)
            varargs = None

            # All args are positional, so convert kargs to positional
            for item in m.annotations:
                annotation = (
                    sig.parameters[item].annotation
                    if item != "return"
                    else sig.return_annotation
                )
                if sys.version_info < (3, 10) and isinstance(annotation, str):
                    annotation = eval(annotation)
                if item == m.varargs:
                    varargs = annotation
                elif item != "return":
                    positional[args_names.index(item)] = annotation

            tailargs = self._positional_validate(
                tailargs, positional, varargs, m.args[1:], m.varargs
            )

        ordered = [
            (f, a)
            for _, f, a in sorted((sf.index, f, sf.val) for f, sf in swfuncs.items())
        ]
        return ordered, tailargs

    def _positional_validate(self, args, validator_list, varargs, argnames, varargname):
        """Makes sure args follows the validation given input"""
        out_args = list(args)

        for i in range(min(len(args), len(validator_list))):
            if validator_list[i] is not None:
                out_args[i] = self._handle_argument(
                    args[i], validator_list[i], argnames[i]
                )

        if len(args) > len(validator_list):
            if varargs is not None:
                out_args[len(validator_list) :] = [
                    self._handle_argument(a, varargs, varargname)
                    for a in args[len(validator_list) :]
                ]
            else:
                out_args[len(validator_list) :] = args[len(validator_list) :]

        return out_args

    @classmethod
    def run(
        cls,
        argv=None,
        exit=True,  # pylint: disable=redefined-builtin
    ):
        """
        Runs the application, taking the arguments from ``sys.argv`` by default if
        nothing is passed. If ``exit`` is
        ``True`` (the default), the function will exit with the appropriate return code;
        otherwise it will return a tuple of ``(inst, retcode)``, where ``inst`` is the
        application instance created internally by this function and ``retcode`` is the
        exit code of the application.

        .. note::
           Setting ``exit`` to ``False`` is intendend for testing/debugging purposes only -- do
           not override it in other situations.
        """
        if argv is None:
            argv = sys.argv
        cls.autocomplete(argv)
        argv = list(argv)
        inst = cls(argv.pop(0))
        retcode = 0
        try:
            swfuncs, tailargs = inst._parse_args(argv)
            ordered, tailargs = inst._validate_args(swfuncs, tailargs)
        except ShowHelp:
            inst.help()
        except ShowHelpAll:
            inst.helpall()
        except ShowVersion:
            inst.version()
        except SwitchError as ex:
            print(T_("Error: {0}").format(ex))
            print(T_("------"))
            inst.help()
            retcode = 2
        else:
            for f, a in ordered:
                f(inst, *a)

            cleanup = None
            if not inst.nested_command or inst.CALL_MAIN_IF_NESTED_COMMAND:
                retcode = inst.main(*tailargs)
                cleanup = functools.partial(inst.cleanup, retcode)
            if not retcode and inst.nested_command:
                subapp, argv = inst.nested_command
                subapp.parent = inst
                inst, retcode = subapp.run(argv, exit=False)

            if cleanup:
                cleanup()

            if retcode is None:
                retcode = 0

        if exit:
            sys.exit(retcode)
        else:
            return inst, retcode

    @classmethod
    def invoke(cls, *args, **switches):
        """Invoke this application programmatically (as a function), in the same way ``run()``
        would. There are two key differences: the return value of ``main()`` is not converted to
        an integer (returned as-is), and exceptions are not swallowed either.

        :param args: any positional arguments for ``main()``
        :param switches: command-line switches are passed as keyword arguments,
                         e.g., ``foo=5`` for ``--foo=5``
        """

        inst = cls("")

        swfuncs = inst._parse_kwd_args(switches)
        ordered, tailargs = inst._validate_args(swfuncs, args)
        for f, a in ordered:
            f(inst, *a)

        cleanup = None
        if not inst.nested_command or inst.CALL_MAIN_IF_NESTED_COMMAND:
            retcode = inst.main(*tailargs)
            cleanup = functools.partial(inst.cleanup, retcode)
        if not retcode and inst.nested_command:
            subapp, argv = inst.nested_command
            subapp.parent = inst
            inst, retcode = subapp.run(argv, exit=False)

        if cleanup:
            cleanup()

        return inst, retcode

    def _parse_kwd_args(self, switches):
        """Parses keywords (positional arguments), used by invoke."""
        swfuncs = {}
        for index, (swname, val) in enumerate(switches.items(), 1):
            switch_local = getattr(type(self), swname)
            swinfo = self._switches_by_func[switch_local._switch_info.func]
            if isinstance(switch_local, CountOf):
                p = (range(val),)
            elif swinfo.list and not hasattr(val, "__iter__"):
                raise SwitchError(
                    T_("Switch {0} must be a sequence (iterable)").format(swname)
                )
            elif not swinfo.argtype:
                # a flag
                if val not in (True, False, None, Flag):
                    raise SwitchError(T_("Switch {0} is a boolean flag").format(swname))
                p = ()
            else:
                p = (val,)
            swfuncs[swinfo.func] = SwitchParseInfo(swname, p, index)
        return swfuncs

    def main(self, *args):
        """Implement me (no need to call super)"""
        if self._subcommands:
            if args:
                print(T_("Unknown sub-command '{0}'").format(args[0]))
                print(T_("------"))
                self.help()
                return 1
            if not self.nested_command:
                print(T_("No sub-command given"))
                print(T_("------"))
                self.help()
                return 1
            return 0

        print(T_("main() not implemented"))
        return 1

    def cleanup(self, retcode):
        """Called after ``main()`` and all sub-applications have executed, to perform any necessary cleanup.

        :param retcode: the return code of ``main()``
        """

    @switch(
        ["--help-all"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints help messages of all sub-commands and quits"""),
    )
    def helpall(self):
        """Prints help messages of all sub-commands and quits"""
        self.help()
        print()

        if self._subcommands:
            for name, subcls in sorted(self._subcommands.items()):
                subapp = (subcls.get())(f"{self.PROGNAME} {name}")
                subapp.parent = self
                for si in subapp._switches_by_func.values():
                    if si.group == "Meta-switches":
                        si.group = "Hidden-switches"
                subapp.helpall()

    @switch(
        ["-h", "--help"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints this help message and quits"""),
    )
    def help(self):  # @ReservedAssignment
        """Prints this help message and quits"""
        if self._get_prog_version():
            self.version()
            print()
        if self.DESCRIPTION:
            print(self.DESCRIPTION.strip() + "\n")

        def split_indentation(s):
            """Identifies the initial indentation (all spaces) of the string and returns the indentation as well
            as the remainder of the line.
            """
            i = 0
            while i < len(s) and s[i] == " ":
                i += 1
            return s[:i], s[i:]

        def paragraphs(text):
            """Yields each paragraph of text along with its initial and subsequent indentations to be used by
            textwrap.TextWrapper.

            Identifies list items from their first non-space character being one of bullets '-', '*', and '/'.
            However, bullet '/' is invisible and is removed from the list item.

            :param text: The text to separate into paragraphs
            """

            paragraph = None
            initial_indent = ""
            subsequent_indent = ""

            def current():
                """Yields the current result if present."""
                if paragraph:
                    yield paragraph, initial_indent, subsequent_indent

            for part in text.lstrip("\n").split("\n"):
                indent, line = split_indentation(part)

                if len(line) == 0:
                    # Starting a new paragraph
                    yield from current()
                    yield "", "", ""

                    paragraph = None
                    initial_indent = ""
                    subsequent_indent = ""
                else:
                    # Adding to current paragraph
                    def is_list_item(line):
                        """Returns true if the first element of 'line' is a bullet character."""
                        bullets = ["-", "*", "/"]
                        return line[0] in bullets

                    def has_invisible_bullet(line):
                        """Returns true if the first element of 'line' is the invisible bullet ('/')."""
                        return line[0] == "/"

                    if is_list_item(line):
                        # Done with current paragraph
                        yield from current()

                        if has_invisible_bullet(line):
                            line = line[1:]

                        paragraph = line
                        initial_indent = indent

                        # Calculate extra indentation for subsequent lines of this list item
                        i = 1
                        while i < len(line) and line[i] == " ":
                            i += 1
                        subsequent_indent = indent + " " * i
                    else:
                        if not paragraph:
                            # Start a new paragraph
                            paragraph = line
                            initial_indent = indent
                            subsequent_indent = indent
                        else:
                            # Add to current paragraph
                            paragraph = paragraph + " " + line

            yield from current()

        def wrapped_paragraphs(text, width):
            """Yields each line of each paragraph of text after wrapping them on 'width' number of columns.

            :param text: The text to yield wrapped lines of
            :param width: The width of the wrapped output
            """
            if not text:
                return

            width = max(width, 1)

            for paragraph, initial_indent, subsequent_indent in paragraphs(text):
                wrapper = TextWrapper(
                    width,
                    initial_indent=initial_indent,
                    subsequent_indent=subsequent_indent,
                )
                w = wrapper.wrap(paragraph)
                yield from w
                if len(w) == 0:
                    yield ""

        cols, _ = get_terminal_size()
        for line in wrapped_paragraphs(self.DESCRIPTION_MORE, cols):
            print(line)

        m = inspect.getfullargspec(self.main)
        tailargs = m.args[1:]  # skip self
        if m.defaults:
            for i, d in enumerate(reversed(m.defaults)):
                tailargs[-i - 1] = f"[{tailargs[-i - 1]}={d}]"
        if m.varargs:
            tailargs.append(f"{m.varargs}...")
        tailargs = " ".join(tailargs)

        utc = self.COLOR_USAGE_TITLE or self.COLOR_USAGE
        print(utc | T_("Usage:"))

        with self.COLOR_USAGE:
            if not self.USAGE:
                if self._subcommands:
                    self.USAGE = T_(
                        "    {progname} [SWITCHES] [SUBCOMMAND [SWITCHES]] {tailargs}\n"
                    )
                else:
                    self.USAGE = T_("    {progname} [SWITCHES] {tailargs}\n")
            print(
                self.USAGE.format(
                    progname=colors.filter(self.PROGNAME), tailargs=tailargs
                )
            )

        by_groups = {}
        for si in self._switches_by_func.values():
            if si.group not in by_groups:
                by_groups[si.group] = []
            by_groups[si.group].append(si)

        def switchs(by_groups, show_groups):
            for grp, swinfos in sorted(by_groups.items(), key=lambda item: item[0]):
                if show_groups:
                    lgrp = T_(grp) if grp in _switch_groups else grp
                    print(self.COLOR_GROUP_TITLES[grp] | lgrp + ":")

                for si in sorted(swinfos, key=lambda x: x.names):
                    swnames = ", ".join(
                        ("-" if len(n) == 1 else "--") + n
                        for n in si.names
                        if n in self._switches_by_name
                        and self._switches_by_name[n] == si
                    )
                    if si.argtype:
                        if hasattr(si.argtype, "__name__"):
                            typename = si.argtype.__name__
                        else:
                            typename = str(si.argtype)
                        argtype = f" {si.argname.upper()}:{typename}"
                    else:
                        argtype = ""
                    prefix = swnames + argtype
                    yield si, prefix, self.COLOR_GROUPS[grp]

                if show_groups:
                    print()

        sw_width = (
            max(len(prefix) for si, prefix, color in switchs(by_groups, False)) + 4
        )
        description_indent = "    {0}{1}{2}"
        wrapper = TextWrapper(width=max(cols - min(sw_width, 60), 50) - 6)
        indentation = "\n" + " " * (cols - wrapper.width)

        for switch_info, prefix, color in switchs(by_groups, True):
            help_txt = switch_info.help
            if switch_info.list:
                help_txt += T_("; may be given multiple times")
            if switch_info.mandatory:
                help_txt += T_("; required")
            if switch_info.requires:
                help_txt += T_("; requires {0}").format(
                    ", ".join(
                        (("-" if len(switch) == 1 else "--") + switch)
                        for switch in switch_info.requires
                    )
                )
            if switch_info.excludes:
                help_txt += T_("; excludes {0}").format(
                    ", ".join(
                        (("-" if len(switch) == 1 else "--") + switch)
                        for switch in switch_info.excludes
                    )
                )

            msg = indentation.join(
                wrapper.wrap(" ".join(ln.strip() for ln in help_txt.splitlines()))
            )

            if len(prefix) + wrapper.width >= cols:
                padding = indentation
            else:
                padding = " " * max(cols - wrapper.width - len(prefix) - 4, 1)
            print(description_indent.format(color | prefix, padding, color | msg))

        if self._subcommands:
            gc = self.COLOR_GROUP_TITLES["Sub-commands"]
            print(gc | T_("Sub-commands:"))
            for name, subcls in sorted(self._subcommands.items()):
                with gc:
                    subapp = subcls.get()
                    doc = subapp.DESCRIPTION or getdoc(subapp)
                    if self.SUBCOMMAND_HELPMSG:
                        help_str = doc + "; " if doc else ""
                        help_str += self.SUBCOMMAND_HELPMSG.format(
                            parent=self.PROGNAME, sub=name
                        )
                    else:
                        help_str = doc or ""

                    msg = indentation.join(
                        wrapper.wrap(
                            " ".join(ln.strip() for ln in help_str.splitlines())
                        )
                    )

                    if len(name) + wrapper.width >= cols:
                        padding = indentation
                    else:
                        padding = " " * max(cols - wrapper.width - len(name) - 4, 1)
                    if colors.contains_colors(subcls.name):
                        bodycolor = colors.extract(subcls.name)
                    else:
                        bodycolor = gc

                    print(
                        description_indent.format(
                            subcls.name, padding, bodycolor | colors.filter(msg)
                        )
                    )

    def _get_prog_version(self):
        ver = None
        curr = self
        while curr is not None:
            ver = getattr(curr, "VERSION", None)
            if ver is not None:
                return ver
            curr = curr.parent
        return ver

    @switch(
        ["-v", "--version"],
        overridable=True,
        group="Meta-switches",
        help=T_("""Prints the program's version and quits"""),
    )
    def version(self):
        """Prints the program's version and quits"""
        ver = self._get_prog_version()
        ver_name = ver if ver is not None else T_("(version not set)")
        print(f"{self.PROGNAME} {ver_name}")
